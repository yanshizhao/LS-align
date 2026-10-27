#include "mini_test.h"
#include "Ligands.h"
#include "LigRotamerDetector.h"
#include "ShortestRoad.h"
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
