#include "mini_test.h"
#include "leak_probe.h"
#include "Ligands.h"
#include "LigRotamerDetector.h"
#include "FlexibleLSalign.h"
#include "RigidLSalign.h"
#include <vector>

/*
 * test_xfail.cpp -- regression locks for defects fixed in the batch-1 round.
 *
 * The file name is kept for compatibility with older references, but every
 * case here is now a normal TEST: the defects it encodes are fixed, so the
 * case must pass.  If one starts failing again, the corresponding defect has
 * regressed.
 *
 * Locked defects:
 *   ligand_operator_assign_does_not_accumulate          -- B9
 *   leak_dijkstra_flag_array                            -- B6
 *   leak_rigid_lsalign                                  -- B1 + B6 + B19
 *   leak_flexible_lsalign_align                         -- B4 + B15 + B17 + B11-B13
 *   aligned_rmsd_nan_when_no_pairs                      -- B8
 */

/* ---- Bug: Ligand::operator= does not clear points() -> accumulates ---- */
TEST(ligand_operator_assign_does_not_accumulate) {
    Ligands ligs("fixtures/self.mol2", false);
    Ligand b;                 /* default ctor */
    b = ligs[0];              /* operator= #1 */
    b = ligs[0];              /* operator= #2 */
    /* correct behaviour: size stays equal to the source */
    CHECK(b.size() == ligs[0].size());
}

/* ---- Bug: DijkstraSRoad leaks the 'flag' array on every getRoad() ---- */
TEST(leak_dijkstra_flag_array) {
    Ligands ligs("fixtures/octane.mol2", false);
    Ligand L = ligs[0];
    /* warm up lazy allocations (locale, iostream) outside the measured scope */
    { std::vector<std::vector<BondType> > warm = L.getAtomBondTypes(); (void)warm; }
    LeakDelta d;
    /* L2 must be inside the measured scope: its bondDisMtx/btMtx/atomBondTypes
     * caches are only freed by ~Ligand, so measuring while L2 is still alive
     * counts those caches as leaks and can never go green on the B6 fix. */
    { Ligand L2 = ligs[0]; L2.getAtomBondTypes(); }
    /* correct behaviour: no net allocation leak */
    CHECK(d.delta() == 0);
}

/* ---- Bug: RigidLSalign destructor never frees p_AliInfo (inverted test) ---- */
TEST(leak_rigid_lsalign) {
    Ligands a("fixtures/self.mol2", false);
    Ligands b("fixtures/copy.mol2", false);
    { LeakDelta warm; { RigidLSalign t(a[0], b[0], 2.0); (void)t; } (void)warm; }
    LeakDelta d;
    { RigidLSalign ali(a[0], b[0], 2.0); (void)ali.getLSscore8QuerySize(); }
    /* correct behaviour: constructing+destroying a RigidLSalign frees everything */
    CHECK(d.delta() == 0);
}

/* ---- Bug: FlexibleLSalign::align() leaks non-winning rotamer objects ---- */
TEST(leak_flexible_lsalign_align) {
    Ligands q("fixtures/octane.mol2", false);
    Ligands t("fixtures/self.mol2", false);
    { LeakDelta warm; { FlexibleLSalign f(q[0], t[0], 2.0, false, false, false); RigidLSalign* p = f.align(10); delete p; } (void)warm; }
    LeakDelta d;
    {
        FlexibleLSalign f(q[0], t[0], 2.0, false, false, false);
        RigidLSalign* p = f.align(10);
        delete p;
    }
    /* correct behaviour: no leak from the rotamer tournament */
    CHECK(d.delta() == 0);
}

/* ---- Bug: getAlignedRMSD() divides by zero when nothing aligns -> NaN ---- */
TEST(aligned_rmsd_nan_when_no_pairs) {
    /* unrelatedA/B are 3 tiny dissimilar atoms; d0 clamps to 0.1 A, so
     * finalAlignedPair is empty and the RMSD computation is 0/0. */
    Ligands a("fixtures/unrelatedA.mol2", false);
    Ligands b("fixtures/unrelatedB.mol2", false);
    RigidLSalign ali(a[0], b[0], 0.1);
    double r = ali.getAlignedRMSD();
    /* correct behaviour: a defined (finite) value, not NaN */
    CHECK(r == r);            /* NaN != NaN */
    CHECK(r >= 0.0);
}
