#include "mini_test.h"
#include "leak_probe.h"
#include "Ligands.h"
#include "RigidLSalign.h"
#include "LigRotamerDetector.h"
#include "ShortestRoad.h"
#include "EnhancedGS.h"
#include "Evaluation.h"
#include "AtomFixInfo.h"
#include <vector>
#include <string>

/*
 * test_delete_mismatch.cpp -- red locks for defects in the allocation /
 * copy-semantics area (B5/B9/B12/B13/B22/B24 + U2).
 *
 * Every XFAIL below encodes the CORRECT behaviour, so it is red before the
 * corresponding fix and turns XPASS afterwards.  The comment above each case
 * states why it is expected to be red today.
 *
 * U2 (default ctor leaves bondDisMtx uninitialized) is NOT a deterministic red
 * lock: its observable effect depends on whatever happens to be on the stack.
 * It is therefore a plain TEST (see the note on that case).
 */

/* ---- B12 + B13: ShortestRoad frees int[1] with scalar delete ---- */
XFAIL(no_delete_mismatch_shortestroad) {
    /* pre-fix: 65 mismatches on octane (64 x B13 in the n^2 loop + 1 x B12) */
    Ligands ligs("fixtures/octane.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];

    int n = L.size();
    double** adj = new2DArr(n, n);
    std::vector<BondInfo> bis = L.getBondInfoes();
    for (std::size_t i = 0; i < bis.size(); ++i) {
        adj[bis[i].aInd][bis[i].bInd] = bis[i].bondDis;
        adj[bis[i].bInd][bis[i].aInd] = bis[i].bondDis;
    }

    MismatchDelta d;
    {
        ShortestRoad sr(adj, n);
        (void)sr.getLongestShortestRoad8Floyd();
    }
    /* correct behaviour: no delete/delete[] mismatch */
    CHECK(d.delta() == 0);

    release2DArr(n, adj);
}

/* ---- B5 + B12 + B13: LigRotamerDetector's ShortestRoad and select_N ---- */
XFAIL(no_delete_mismatch_rotamer_c13) {
    /* pre-fix: 170 mismatches on the C13 chain (169 in the detector ctor --
     * ShortestRoad, plus 1 x B5 inside generateRotamers()).  Reading delta()
     * before the detector is destroyed skips its own B12 (+1). */
    Ligands ligs("fixtures/chain_c13.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];

    MismatchDelta d;
    LigRotamerDetector rd(L);
    std::vector<Ligand> rot = rd.generateRotamers();
    (void)rot;

    /* correct behaviour: no delete/delete[] mismatch */
    CHECK(d.delta() == 0);
}

/* ---- B11: EnhancedGS deletes each int[] with scalar delete ---- */
XFAIL(no_delete_mismatch_enhancedgs) {
    /* pre-fix: one mismatch per row of the score matrix (3 for a 3x3) */
    std::vector<std::vector<double> > m;
    for (int i = 0; i < 3; ++i) {
        std::vector<double> row(3, 0.0);
        row[i] = 1.0;
        m.push_back(row);
    }

    MismatchDelta d;
    {
        EnhancedGS egs(m);
        (void)egs.getPAliInfo()->sco;
    }
    /* correct behaviour: no delete/delete[] mismatch */
    CHECK(d.delta() == 0);
}

/* ---- B22: calculateThreeSpecicalDistPow2 reads [-1] when degenerate ---- */
XFAIL(three_special_dists_degenerate) {
    /* pre-fix: lig_size=1 leaves AInd=BInd=-1, so the second loop reads
     * p_disPow2Mtx[0][-1].  With a plain new2DArr(1,1) that is adjacent heap
     * memory and the result is flaky (sometimes 0, sometimes ~1e-307).  We
     * therefore place a known guard value in the word just before row 0: the
     * pre-fix code deterministically returns {0, 1, 1}; the fixed code never
     * reads the matrix for lig_size < 2 and returns {0, 0, 0}. */
    double* guard_and_cell = new double[2];
    guard_and_cell[0] = 1.0;            /* read as m[0][-1] before the fix */
    guard_and_cell[1] = 0.0;            /* logical m[0][0] */
    double* row = guard_and_cell + 1;
    double** m1 = &row;

    double* r1 = calculateThreeSpecicalDistPow2(1, m1);
    /* correct behaviour: a single atom has no distances -> all zero */
    CHECK(r1[0] == 0.0);
    CHECK(r1[1] == 0.0);
    CHECK(r1[2] == 0.0);
    delete[] r1;
    delete[] guard_and_cell;

    /* two distinct atoms: the A-B distance is the only non-zero reference */
    double** m2 = new2DArr(2, 2);
    m2[0][1] = 2.25;
    m2[1][0] = 2.25;
    double* r2 = calculateThreeSpecicalDistPow2(2, m2);
    CHECK(r2[0] == 2.25);
    CHECK(r2[1] == 0.0);
    CHECK(r2[2] == 0.0);
    delete[] r2;
    release2DArr(2, m2);
}

/* ---- B24: RMSD_lb has no sentinel for an empty input ---- */
XFAIL(rmsd_lb_empty_returns_sentinel) {
    /* pre-fix: (empty, one) = 316227766.0..., (one, empty) = NaN */
    std::vector<Point3D> empty;
    std::vector<Point3D> one;
    one.push_back(Point3D(0.0, 0.0, 0.0));

    /* correct behaviour: a defined sentinel for "not computable" */
    CHECK(RMSD_lb(empty, one) == -1.0);
    CHECK(RMSD_lb(one, empty) == -1.0);
}

/* ---- U2: default ctor does not initialize bondDisMtx ----
 * U2: latent - default-constructed bondDisMtx is garbage but usually NULL on this stack; not deterministically red pre-fix
 * Plain TEST (not XFAIL): the bug is only observable through the destructor and
 * depends on stack contents, so a red lock here would be flaky.  This case
 * records the intended contract and becomes a regression lock once T-015
 * initializes both matrices.
 *
 * NOTE: a bare "{ Ligand b; }" is NOT safe to write here -- pre-fix its
 * destructor deterministically reads the garbage bondDisMtx and the process
 * aborts (exit 127, 3/3 runs), which would kill the whole suite.  The
 * copy/assign sequence below overwrites the garbage pointer with a valid
 * zero-length array, so b can be destroyed at end of scope deterministically
 * (5/5 runs) both before and after the fix.
 */
TEST(ligand_default_ctor_safe_to_destroy) {
    Ligand b;
    Ligand c = b;
    b = c;
    CHECK(b.size() == 0);
    /* b (and c) are destroyed here, at the end of the case's scope */
}

/* ---- B9: operator= appends points instead of replacing them ---- */
XFAIL(ligand_assign_no_accumulate) {
    /* pre-fix: two assignments into a default target leave 4 atoms, not 2 */
    std::vector<Point3D> pts;
    pts.push_back(Point3D(0.0, 0.0, 0.0));
    pts.push_back(Point3D(1.5, 0.0, 0.0));
    std::vector<std::string> ats;
    ats.push_back("C");
    ats.push_back("C");
    std::vector<BondInfo> bis;
    BondInfo bi;
    bi.aInd = 0;
    bi.bInd = 1;
    bi.bt = SINGLE;
    bi.bondDisPow2 = 2.25;
    bi.bondDis = 1.5;
    bis.push_back(bi);
    Ligand a("TWO", pts, ats, bis);

    Ligand c;
    c = a;
    c = a;

    /* correct behaviour: assignment replaces the target's contents */
    CHECK(c.size() == a.size());
}

/* ---- B9: operator= leaves stale atomBondTypes in a non-empty target ---- */
XFAIL(ligand_copy_rebuilds_identical_atom_bond_types) {
    /* pre-fix: dst keeps octane's 8-atom cache and octane's 8 points, so
     * dst.size()==14 and atomBondTypes.size()==8 after dst = self[0]. */
    Ligands octLigs("fixtures/octane.mol2", false);
    Ligands selfLigs("fixtures/self.mol2", false);
    std::vector<Ligand> oct = octLigs.getLigands();
    std::vector<Ligand> self = selfLigs.getLigands();

    Ligand dst = oct[0];                 /* 8 atoms */
    dst.getAtomBondTypes();              /* make dst's cache non-empty */
    Ligand src = self[0];                /* 6 atoms */

    dst = src;

    /* correct behaviour: dst is a faithful copy of src, caches rebuilt */
    CHECK(dst.size() == 6);
    CHECK(dst.getAtomBondTypes().size() == 6);
}

/* ---- B18: Ra/Re entries are swapped/duplicated in the mass and vdW tables ---- */
XFAIL(mass_vdw_table_ra_re) {
    /* Pre-fix: mass table has Ra=186.2071 / Re=226 (swapped), vdW has Re=2.43
     * (duplicate of the correct Re=2.05 earlier in the string).  Expected
     * values come from the periodic table: Ra 226.0254, Re 186.207. */
    AtomRelativeMass mass;
    AtomVarDerWaalRadius vdw;
    CHECK_CLOSE(mass[std::string("Ra")], 226.0, 1.0);
    CHECK_CLOSE(mass[std::string("Re")], 186.2071, 0.01);
    CHECK_CLOSE(vdw[std::string("Re")], 2.05, 0.01);
    CHECK_CLOSE(vdw[std::string("Ra")], 2.43, 0.01);
}
