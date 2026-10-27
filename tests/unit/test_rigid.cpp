#include "mini_test.h"
#include "Ligands.h"
#include "RigidLSalign.h"

/* Rigid-LS-align: self-alignment and rigid-transform invariance.
 * For a perfect alignment every pair contributes 0.45 + 0.10 + 0.45 = 1.0,
 * so getLSscore8QuerySize() == sum/Nq -> ~1.0. */

TEST(rigid_self_alignment_scores_one) {
    Ligands ligs("fixtures/self.mol2", false);
    RigidLSalign ali(ligs[0], ligs[0], 2.0);
    CHECK_CLOSE(ali.getLSscore8QuerySize(), 1.0, 1e-3);
    CHECK_CLOSE(ali.getLSscore8TemplSize(), 1.0, 1e-3);
    CHECK_CLOSE(ali.getJaccardRatio(), 1.0, 1e-3);
    CHECK_CLOSE(ali.getRMSD_lb(), 0.0, 1e-3);
}

TEST(rigid_transform_invariance) {
    /* copy.mol2 is self.mol2 rigidly rotated+translated -> still PC ~ 1.0 */
    Ligands a("fixtures/self.mol2", false);
    Ligands b("fixtures/copy.mol2", false);
    RigidLSalign ali(a[0], b[0], 2.0);
    CHECK_CLOSE(ali.getLSscore8QuerySize(), 1.0, 1e-3);
    CHECK_CLOSE(ali.getJaccardRatio(), 1.0, 1e-3);
    CHECK_CLOSE(ali.getRMSD_lb(), 0.0, 1e-3);
}

TEST(rigid_ace_example_matches_paper_band) {
    /* Paper (Sec 3.1): ZINC03812885 vs the ACE seed 3bkl.
     * The paper reports an alignment of ~13 pairs, RMSD ~= 0.44 A, LS-score ~= 0.62.
     * We assert loose bands because the shipped CLI metrics differ (RMSD_lb is a
     * lower-bound RMSD, not the aligned RMSD; PC-score8Q is size-normalised). */
    Ligands q("../example/ace_ligands.mol2", false);
    Ligands t("../example/xtal-lig.mol2", false);
    CHECK(q.size() == 49);
    CHECK(t.size() == 1);

    /* locate ZINC03812885 among the 49 queries */
    int idx = -1;
    for (int i = 0; i < q.size(); ++i) {
        if (q[i].getLigname() == "ZINC03812885") { idx = i; break; }
    }
    CHECK(idx >= 0);
    if (idx < 0) return;

    RigidLSalign ali(q[idx], t[0], 2.0);
    double lsScore = ali.getLSscore8TemplSize();  /* vs 29-atom template */
    /* loose band around the paper's ~0.62 */
    CHECK(lsScore > 0.30 && lsScore < 0.95);
    /* it should align a healthy number of atom pairs */
    CHECK(ali.getFinalAlignedPair().size() >= 5);
}
