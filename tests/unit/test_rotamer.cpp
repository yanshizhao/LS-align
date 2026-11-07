#include "mini_test.h"
#include "Ligands.h"
#include "ShortestRoad.h"
#include "RigidLSalign.h"
#define private public
#include "LigRotamerDetector.h"
#undef private
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

/* LigRotamerDetector: conformer generation by rotating rotatable single bonds.
 * Uses a fixed small fixture (octane) so the deterministic path is exercised
 * (N=3 and at most 11 angles -> at most 1331 conformers, below the 5000 cap). */

/* R5 removed the no-arg generateRotamers() (N=5, no aligned-bond filter).
 * The cases below used it only as "generate the conformers of this ligand";
 * an initial rigid alignment with no aligned pairs reproduces that exactly
 * (every single bond stays in the rotation set). */
static RigidLSalign emptyRigidAlignment() {
    Ligands a("fixtures/tiny.mol2", false);
    Ligands b("fixtures/unrelatedB.mol2", false);
    std::vector<Ligand> av = a.getLigands();
    std::vector<Ligand> bv = b.getLigands();
    return RigidLSalign(av[0], bv[0], 0.1);
}

TEST(rotamer_generation_nonempty_and_size_preserved) {
    Ligands ligs("fixtures/octane.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];
    CHECK(L.size() == 8);

    RigidLSalign rigid = emptyRigidAlignment();
    CHECK(rigid.getFinalAlignedPair().empty());

    LigRotamerDetector rd(L);
    std::vector<Ligand> rot = rd.generateRotamers(&rigid);

    CHECK(rot.size() >= 1);
    for (std::size_t i = 0; i < rot.size(); ++i) {
        /* rotating bonds never changes the atom count */
        CHECK(rot[i].size() == L.size());
        CHECK(rot[i].getLigname() == L.getLigname());
    }
}

TEST(shortestroad_longest_path_on_octane) {
    /* The longest shortest-path in a linear chain of 8 atoms has 8 nodes. */
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
    ShortestRoad sr(adj, n);
    std::vector<int> road = sr.getLongestShortestRoad8Floyd();
    CHECK(road.size() == 8);
    release2DArr(n, adj);
}

/* B26: the uniform-distribution sampler used to read
 * sortedSIInds[i*step + i] unconditionally, which is out of range for
 * i >= mod when i*step + i >= s (chain_c13: s=6; i=3 reads [6], i=4 reads
 * [8] for N=5).  The out-of-range value was always overwritten by the
 * i*step + mod pick before use, so the selected indices are unchanged.
 * This helper exposes the private sampler; it must never read out of range. */
static std::vector<int> b26SelectC13(int N, int& bondCount) {
    Ligands ligs("fixtures/chain_c13.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    LigRotamerDetector rd(all[0]);
    std::vector<BondInfo> singleBonds = rd.findAllSingleBond();
    std::vector<std::vector<std::vector<int> > > twoPartsIndss =
        rd.getUsefulBondTwoPartIndsInfo(&singleBonds);
    bondCount = (int)singleBonds.size();
    CHECK((int)twoPartsIndss.size() == bondCount);
    return rd.select_N_RotableSingleBondOnUniformDistribution(
        twoPartsIndss, singleBonds, N);
}

TEST(b26_c13_no_abort) {
    /* With -D_GLIBCXX_ASSERTIONS the pre-fix read aborts inside
     * vector::operator[]; the fixed sampler must return N in-range indices. */
    int s = 0;
    std::vector<int> sel3 = b26SelectC13(3, s);
    CHECK(s == 6);
    CHECK(sel3.size() == 3);
    for (std::size_t i = 0; i < sel3.size(); ++i)
        CHECK(sel3[i] >= 0 && sel3[i] < s);

    std::vector<int> sel5 = b26SelectC13(NRS_LEGACY_UNPARAM_BOND_NUM, s);
    CHECK(s == 6);
    CHECK(sel5.size() == NRS_LEGACY_UNPARAM_BOND_NUM);
    for (std::size_t i = 0; i < sel5.size(); ++i)
        CHECK(sel5[i] >= 0 && sel5[i] < s);
}

TEST(b26_indices_identical) {
    /* Values measured on chain_c13 (s=6): pre-fix and post-fix sequences are
     * bit-identical because the out-of-range read was always discarded. */
    int s = 0;
    std::vector<int> sel3 = b26SelectC13(3, s);
    CHECK(sel3.size() == 3);
    CHECK(sel3.size() == 3 && sel3[0] == 3 && sel3[1] == 2 && sel3[2] == 0);

    std::vector<int> sel5 = b26SelectC13(5, s);
    CHECK(sel5.size() == 5);
    CHECK(sel5.size() == 5 && sel5[0] == 3 && sel5[1] == 2
          && sel5[2] == 1 && sel5[3] == 0 && sel5[4] == 5);
}

/* F8: the magic numbers are named; the values must stay identical to the
 * literals they replaced (behaviour-preserving refactor). */
TEST(f8_macro_values) {
    CHECK(NRS_ROTATABLE_BOND_NUM == 3);
    CHECK(NRS_MIN_ROTATABLE_BOND_NUM == 3);
    CHECK(NRS_LEGACY_UNPARAM_BOND_NUM == 5);
    CHECK(ROTATABLE_BOND_MIN_SIDE_ATOMS == 4);
}

/* F8: every index the sampler returns must be a valid index into singleBonds,
 * for both the N=3 (main path) and N=5 (legacy unparameterised value; after
 * R5 only these sampler tests still pass it) calls. */
TEST(f8_chain_c13_index_range) {
    int s = 0;
    std::vector<int> sel3 = b26SelectC13(3, s);
    CHECK(s == 6);
    CHECK(sel3.size() == 3);
    for (std::size_t i = 0; i < sel3.size(); ++i)
        CHECK(sel3[i] >= 0 && sel3[i] < s);

    std::vector<int> sel5 = b26SelectC13(5, s);
    CHECK(sel5.size() == 5);
    for (std::size_t i = 0; i < sel5.size(); ++i)
        CHECK(sel5[i] >= 0 && sel5[i] < s);

    /* the sampler clamps N up to the minimum (here the ligand yields s=6) */
    std::vector<int> sel0 = b26SelectC13(0, s);
    CHECK(sel0.size() == NRS_MIN_ROTATABLE_BOND_NUM);
    for (std::size_t i = 0; i < sel0.size(); ++i)
        CHECK(sel0[i] >= 0 && sel0[i] < s);
}

/* ================= B31 + B33 (ShortestRoad robustness) ================= */

static double** b33AdjFromLigand(Ligand& L) {
    int n = L.size();
    double** adj = new2DArr(n, n);
    std::vector<BondInfo> bis = L.getBondInfoes();
    for (std::size_t i = 0; i < bis.size(); ++i) {
        adj[bis[i].aInd][bis[i].bInd] = bis[i].bondDis;
        adj[bis[i].bInd][bis[i].aInd] = bis[i].bondDis;
    }
    return adj;
}

/* A simple in-memory ligand: n C atoms and the given single bonds. */
static Ligand b33MakeLigand(int n, const std::vector<BondInfo>& bonds) {
    std::vector<Point3D> pts;
    std::vector<std::string> types;
    for (int i = 0; i < n; ++i) {
        pts.push_back(Point3D(i * 1.5, (i % 2) * 1.1, 0.0));
        types.push_back("C");
    }
    return Ligand("FRAGS", pts, types, bonds);
}

/* B33: an unreachable pair must not be encoded as a direct edge.  Pre-fix the
 * path was the fake two-element [0, 3]; it now holds only the start atom. */
TEST(b33_hop_unreachable_isolated) {
    Ligands ligs("fixtures/isolated_atom.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];
    int n = L.size();
    CHECK(n == 4);

    double** adj = b33AdjFromLigand(L);
    ShortestRoad sr(adj, n);
    CHECK(sr.isValid());
    CHECK(sr.getHopCount(0, 1) == 1);    /* direct bond */
    CHECK(sr.getHopCount(0, 2) == 2);    /* two bonds */
    CHECK(sr.getHopCount(0, 3) == -1);   /* isolated atom: unreachable */
    CHECK(sr.getHopCount(3, 0) == -1);
    CHECK(sr.getHopCount(1, 1) == 0);

    int len = 0;
    while (len < n && sr.getAllRoad8Floyd()[0][3][len] >= 0) ++len;
    CHECK(len == 1);                     /* only the start atom, not [0, 3] */

    /* out-of-range indices are rejected instead of dereferenced */
    CHECK(sr.getHopCount(-1, 0) == -1);
    CHECK(sr.getHopCount(0, n) == -1);
    release2DArr(n, adj);
}

/* B33: connected molecules keep the exact hop counts. */
TEST(b33_hop_octane_connected) {
    Ligands ligs("fixtures/octane.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];
    int n = L.size();
    CHECK(n == 8);

    double** adj = b33AdjFromLigand(L);
    ShortestRoad sr(adj, n);
    CHECK(sr.isValid());
    CHECK(sr.getHopCount(0, 1) == 1);
    CHECK(sr.getHopCount(0, 7) == 7);
    CHECK(sr.getHopCount(7, 0) == 7);
    CHECK(sr.getHopCount(4, 4) == 0);
    for (int i = 0; i + 1 < n; ++i)
        CHECK(sr.getHopCount(i, i + 1) == 1);
    release2DArr(n, adj);
}

/* B31: the n^3 budget is checked before any allocation.  A rejected build must
 * leave every pointer NULL and allocate nothing at all (108 MB at n=301). */
TEST(b31_shortestroad_rejects_oversized_without_alloc) {
    CHECK(ShortestRoad::MAX_FLOYD_NODES == 300);

    double** adj = new2DArr(2, 2);
    adj[0][1] = adj[1][0] = 1.5;

    {
        LeakDelta d;
        ShortestRoad big(adj, 301);
        CHECK(d.delta() == 0);           /* rejected before new2DArr/new3DIntArr */
        CHECK(!big.isValid());
        CHECK(big.getNodeNum() == 0);
        CHECK(big.getHopCount(0, 0) == -1);
        CHECK(big.getLongestShortestRoad8Floyd().empty());
        CHECK(big.getAllRoad8Floyd() == (int***)NULL);
    }
    {
        ShortestRoad zero(adj, 0);
        CHECK(!zero.isValid());
        CHECK(zero.getNodeNum() == 0);
    }
    {
        ShortestRoad neg(adj, -3);
        CHECK(!neg.isValid());
        CHECK(neg.getHopCount(0, 0) == -1);
    }
    release2DArr(2, adj);
}

/* B31: the detector skips an oversized ligand (no bonds returned) and prints
 * one diagnostic; it must not crash or allocate the Floyd tables. */
TEST(b31_oversized_ligand_skipped_with_diagnostic) {
    BondInfo b;
    b.aInd = 0;
    b.bInd = 1;
    b.bt = SINGLE;
    b.bondDis = 1.5;
    b.bondDisPow2 = 2.25;
    std::vector<BondInfo> bonds;
    bonds.push_back(b);
    Ligand L = b33MakeLigand(301, bonds);
    CHECK(L.size() == 301);

    RigidLSalign rigid = emptyRigidAlignment();
    std::ostringstream captured;
    std::streambuf* oldBuf = std::cerr.rdbuf(captured.rdbuf());
    LigRotamerDetector rd(L);
    std::vector<Ligand> rot = rd.generateRotamers(&rigid);
    std::cerr.rdbuf(oldBuf);

    CHECK(rot.empty());
    CHECK(captured.str().find("above the supported maximum") != std::string::npos);
}

/* B33 (both layers): 9-bond chain plus one fully isolated atom.  Every bond is
 * rejected (no published partition), no rotamers are generated and the old
 * data-loss diagnostic ERROR: 201702161138 is never printed. */
TEST(b33_long_chain_isolated_no_rotamers_no_error) {
    Ligands ligs("fixtures/long_chain_isolated.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];
    CHECK(L.size() == 11);

    double** adj = b33AdjFromLigand(L);
    {
        ShortestRoad sr(adj, L.size());
        CHECK(sr.getHopCount(0, 9) == 9);
        CHECK(sr.getHopCount(0, 10) == -1);
        CHECK(sr.getHopCount(10, 9) == -1);
    }
    release2DArr(L.size(), adj);

    RigidLSalign rigid = emptyRigidAlignment();
    std::ostringstream captured;
    std::streambuf* oldBuf = std::cerr.rdbuf(captured.rdbuf());
    LigRotamerDetector rd(L);
    std::vector<BondInfo> singleBonds = rd.findAllSingleBond();
    bool allPartitionsEmpty = true;
    for (std::size_t k = 0; k < singleBonds.size(); ++k) {
        std::vector<std::vector<int> > parts = rd.leftRightPartOnBond(singleBonds[k]);
        if (!parts.empty()) allPartitionsEmpty = false;
    }
    LigRotamerDetector rd2(L);
    std::vector<Ligand> rot = rd2.generateRotamers(&rigid);
    std::cerr.rdbuf(oldBuf);

    CHECK(singleBonds.size() == 9);
    CHECK(allPartitionsEmpty);
    CHECK(rot.empty());
    CHECK(captured.str().find("ERROR: 201702161138") == std::string::npos);
}

/* B33 layer 2: the two sides must cover every atom.  Seven atoms but only two
 * disjoint C-C bonds (0-3 and 1-2); for bond 0-3 the split is {0} vs {3} with
 * atoms {1,2,4,5,6} in neither group.  The old code published that 2/7 split;
 * the coverage check now refuses it at the source. */
TEST(b33_partition_refuses_incomplete_cover) {
    BondInfo b0;
    b0.aInd = 0; b0.bInd = 3; b0.bt = SINGLE;
    b0.bondDis = 4.5; b0.bondDisPow2 = 20.25;
    BondInfo b1;
    b1.aInd = 1; b1.bInd = 2; b1.bt = SINGLE;
    b1.bondDis = 2.2; b1.bondDisPow2 = 4.84;
    std::vector<BondInfo> bonds;
    bonds.push_back(b0);
    bonds.push_back(b1);
    Ligand L = b33MakeLigand(7, bonds);

    LigRotamerDetector rd(L);
    std::vector<BondInfo> singleBonds = rd.findAllSingleBond();
    CHECK(singleBonds.size() == 2);
    for (std::size_t k = 0; k < singleBonds.size(); ++k) {
        std::vector<std::vector<int> > parts = rd.leftRightPartOnBond(singleBonds[k]);
        CHECK(parts.empty());
    }
    RigidLSalign rigid = emptyRigidAlignment();
    CHECK(rd.generateRotamers(&rigid).empty());
}

/* ================= R4 (aligned-bond filter) ================= */

/* A ligand built from atoms [first..last] of another ligand, with the bonds
 * whose two ends both fall in that window (indices rebased to 0).  Used to make
 * a partial template so that only part of the query is rigidly aligned. */
static Ligand r4SubLigand(Ligand& src, int first, int last) {
    std::vector<Point3D> pts;
    std::vector<std::string> types;
    for (int i = first; i <= last; ++i) {
        pts.push_back(src.getIthPoints(i));
        types.push_back(src.getAtomType(i));
    }
    std::vector<BondInfo> all = src.getBondInfoes();
    std::vector<BondInfo> bonds;
    for (std::size_t i = 0; i < all.size(); ++i) {
        if (all[i].aInd >= first && all[i].aInd <= last
            && all[i].bInd >= first && all[i].bInd <= last) {
            BondInfo b = all[i];
            b.aInd -= first;
            b.bInd -= first;
            bonds.push_back(b);
        }
    }
    return Ligand("R4SUB", pts, types, bonds);
}

/* R4: the rigid-aligned-bond filter used to assign isAIndAligned in the B-atom
 * branch, so isBIndAligned stayed false and no bond was ever removed.  With
 * chain13 aligned onto its first 9 atoms, bonds 0-1 .. 7-8 have both ends
 * aligned at distance 0 and must be dropped; only bonds with an unaligned end
 * survive.  Measured: bug build -> singleBonds 12, 26 rotamers which rotate
 * the aligned core; fixed build -> singleBonds 4, 2 rotamers that move only
 * atoms 9-12, leaving every aligned atom (0..8) in place. */
TEST(r4_aligned_bonds_are_removed) {
    Ligands ligs("fixtures/chain_c13.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand q = all[0];
    CHECK(q.size() == 13);
    Ligand t = r4SubLigand(q, 0, 8);

    RigidLSalign rigid(q, t, 2.0);
    std::vector<AlignedPair> pairs = rigid.getFinalAlignedPair();
    CHECK(pairs.size() == 9);                    /* atoms 0..8, distance 0 */

    LigRotamerDetector rd(q);
    CHECK(rd.findAllSingleBond().size() == 12);

    std::vector<Ligand> rot = rd.generateRotamers(&rigid);
    CHECK(rot.size() == 2);                      /* measured, not estimated */
    for (std::size_t r = 0; r < rot.size(); ++r) {
        for (int k = 0; k <= 8; ++k) {
            /* the aligned core must not be rotated any more */
            CHECK(q.getIthPoints(k).dist(rot[r].getIthPoints(k)) < 1e-9);
        }
    }
}

/* R4: a bond with only one end aligned stays in the rotation set.  chain_c13
 * is a straight chain, so any atom prefix with the bonds that fit is another
 * chain.  Query = 10-atom chain, template = its first 5 atoms: atoms 0..4 are aligned at distance 0, so bonds 0-1..3-4 are dropped
 * but bond 4-5 (only atom 4 aligned) must survive.  Measured: fixed build keeps
 * 5 bonds (4-5..8-9) and yields 8 rotamers; some of those rotate the side of
 * bond 4-5, i.e. atoms 0..3 move.  No other kept bond (5-6..8-9) can move an
 * atom below 5, so any motion there proves the half-aligned bond was kept.
 * The buggy build keeps all 9 bonds and yields 26 rotamers. */
TEST(r4_half_aligned_bonds_are_kept) {
    Ligands ligs("fixtures/chain_c13.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand src = all[0];
    CHECK(src.size() == 13);

    std::vector<Point3D> pts;
    std::vector<std::string> types;
    for (int i = 0; i < 10; ++i) {
        pts.push_back(src.getIthPoints(i));
        types.push_back(src.getAtomType(i));
    }
    std::vector<BondInfo> srcBonds = src.getBondInfoes();
    std::vector<BondInfo> bonds;
    for (std::size_t i = 0; i < srcBonds.size(); ++i)
        if (srcBonds[i].aInd < 10 && srcBonds[i].bInd < 10)
            bonds.push_back(srcBonds[i]);
    Ligand q("R4C10", pts, types, bonds);
    CHECK(q.size() == 10);
    Ligand t = r4SubLigand(q, 0, 4);

    RigidLSalign rigid(q, t, 2.0);
    CHECK(rigid.getFinalAlignedPair().size() == 5);
    LigRotamerDetector rd(q);
    CHECK(rd.findAllSingleBond().size() == 9);

    std::vector<Ligand> rot = rd.generateRotamers(&rigid);
    CHECK(rot.size() == 8);                      /* measured, not estimated */

    bool lowAtomMoved = false;
    for (std::size_t r = 0; r < rot.size(); ++r) {
        for (int k = 0; k <= 3; ++k) {
            if (q.getIthPoints(k).dist(rot[r].getIthPoints(k)) > 1e-9)
                lowAtomMoved = true;
        }
    }
    CHECK(lowAtomMoved);                         /* bond 4-5 survived */
}

/* ================= X1 (deterministic decimation above 5000) ================= */

static bool x1SameRotamerSets(const std::vector<Ligand>& a, const std::vector<Ligand>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].size() != b[i].size()) return false;
        for (int k = 0; k < a[i].size(); ++k) {
            Point3D pa = a[i].getIthPoints(k);
            Point3D pb = b[i].getIthPoints(k);
            if (pa[0] != pb[0] || pa[1] != pb[1] || pa[2] != pb[2]) return false;
        }
    }
    return true;
}

/* X1: the old subsetting used srand(clock())+rand(), so repeated runs differed.
 * R5 deleted the no-arg path that carried the >5000 case; the with-arg path
 * keeps the same deterministic decimation, but with N=3 and at most 11 detail
 * angles it can never exceed 11^3 = 1331 conformers, so the >5000 branch is
 * UNREACHABLE from any public caller (and from the CLI).  This case therefore
 * locks (a) reproducibility of the with-arg path and (b) the decimation
 * arithmetic itself, driving the branch through a white-box extended angle set
 * (11 + 9 = 20 angles -> 20^3 = 8000 > 5000).  The extension is a test-only
 * instrument; it does not claim the branch is publicly reachable. */
TEST(x1_deterministic_decimation) {
    Ligands ligs("fixtures/chain_c13.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand src = all[0];
    CHECK(src.size() == 13);

    std::vector<Point3D> pts;
    std::vector<std::string> types;
    for (int i = 0; i < 11; ++i) {
        pts.push_back(src.getIthPoints(i));
        types.push_back(src.getAtomType(i));
    }
    std::vector<BondInfo> allBonds = src.getBondInfoes();
    std::vector<BondInfo> bonds;
    for (std::size_t i = 0; i < allBonds.size() && (int)i < 10; ++i)
        bonds.push_back(allBonds[i]);
    Ligand L("C11CHAIN", pts, types, bonds);
    CHECK(L.size() == 11);

    RigidLSalign rigid = emptyRigidAlignment();

    /* (a) public path: 11 detail angles, N=3 -> 11^3 = 1331, no decimation.
     * Three runs must be identical.  Measured: 342 non-overlapping conformers
     * (the rest are rejected by judgeIsOverlapAtomsInLig). */
    LigRotamerDetector rd1(L, true, true);
    std::vector<Ligand> r1 = rd1.generateRotamers(&rigid);
    LigRotamerDetector rd2(L, true, true);
    std::vector<Ligand> r2 = rd2.generateRotamers(&rigid);
    LigRotamerDetector rd3(L, true, true);
    std::vector<Ligand> r3 = rd3.generateRotamers(&rigid);

    CHECK(!r1.empty());
    CHECK(x1SameRotamerSets(r1, r2));
    CHECK(x1SameRotamerSets(r1, r3));
    CHECK(r1.size() == 342);            /* measured, not estimated */
    CHECK(r1.size() <= 5000);           /* below the decimation cap */

    /* (b) white-box: 20 angles -> 8000 > 5000, so the decimation branch runs.
     * Two runs must still be identical.  Measured: 1139 conformers. */
    LigRotamerDetector d1(L, true, true);
    for (int e = 0; e < 9; ++e)
        d1.rotAngles.push_back((25.0 + 20.0 * e) / 360.0 * DOUBLE_PI);
    CHECK(d1.rotAngles.size() == 20);
    std::vector<Ligand> q1 = d1.generateRotamers(&rigid);

    LigRotamerDetector d2(L, true, true);
    for (int e = 0; e < 9; ++e)
        d2.rotAngles.push_back((25.0 + 20.0 * e) / 360.0 * DOUBLE_PI);
    std::vector<Ligand> q2 = d2.generateRotamers(&rigid);

    CHECK(!q1.empty());
    CHECK(x1SameRotamerSets(q1, q2));
    CHECK(q1.size() == 1139);           /* measured, not estimated */
    CHECK(q1.size() <= 5000);           /* decimation cap holds */

    /* the decimation arithmetic for 8000 records capped at 5000: iter 1..4999
     * maps to floor(iter * 8000 / 5000); the map is strictly increasing (the
     * step is 1.6 records), starts at 1 and ends at 7998. */
    {
        const int allRot = 8000;
        const int maxIter = 5000;
        int prev = -1;
        bool strictlyIncreasing = true;
        bool inRange = true;
        for (int iter = 1; iter < maxIter; ++iter) {
            int recordNum = (int)((double)iter * allRot / maxIter);
            if (recordNum >= allRot) recordNum = allRot - 1;
            if (recordNum <= 0 || recordNum >= allRot) inRange = false;
            if (recordNum <= prev) strictlyIncreasing = false;
            prev = recordNum;
        }
        CHECK(inRange);
        CHECK(strictlyIncreasing);
        int firstRec = (int)((double)1 * allRot / maxIter);
        int lastRec = (int)((double)(maxIter - 1) * allRot / maxIter);
        CHECK(firstRec == 1);
        CHECK(lastRec == 7998);
    }
}
