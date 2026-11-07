#include "mini_test.h"
#include "Ligands.h"
#include "ShortestRoad.h"
#include "RigidLSalign.h"
#define private public
#include "LigRotamerDetector.h"
#undef private
#include <vector>

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
