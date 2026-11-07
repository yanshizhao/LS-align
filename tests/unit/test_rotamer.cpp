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
 * (rotamer count 5^N < 5000 -> no srand()-driven random subsetting). */

TEST(rotamer_generation_nonempty_and_size_preserved) {
    Ligands ligs("fixtures/octane.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];
    CHECK(L.size() == 8);

    LigRotamerDetector rd(L);
    std::vector<Ligand> rot = rd.generateRotamers();

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

    std::vector<int> sel5 = b26SelectC13(5, s);
    CHECK(s == 6);
    CHECK(sel5.size() == 5);
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
 * for both the N=3 (main path) and N=5 (legacy unparameterised path) calls. */
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

    std::ostringstream captured;
    std::streambuf* oldBuf = std::cerr.rdbuf(captured.rdbuf());
    LigRotamerDetector rd(L);
    std::vector<Ligand> rot = rd.generateRotamers();
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
    std::vector<Ligand> rot = rd2.generateRotamers();
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
    CHECK(rd.generateRotamers().empty());
}
