#include "mini_test.h"
#include "EnhancedGS.h"
#include <vector>

/* EnhancedGS is the "enhanced greedy search" assignment solver.
 * It maximizes the sum of selected entries with each row/col used at most once.
 * We hand-build matrices whose optimum is obvious. */

TEST(enhancedgs_single_max_is_selected) {
    std::vector<std::vector<double> > m;
    std::vector<double> r0; r0.push_back(0.1); r0.push_back(0.9);
    std::vector<double> r1; r1.push_back(0.8); r1.push_back(0.2);
    m.push_back(r0); m.push_back(r1);

    EnhancedGS egs(m);
    AlignInfo* p = egs.getPAliInfo();
    std::vector<AlignedPair> pairs = p->getAliPairs();
    /* optimum: (0,1)=0.9 and (1,0)=0.8 */
    CHECK(pairs.size() == 2);
    CHECK_CLOSE(p->sco, 1.7, 1e-9);
}

TEST(enhancedgs_identity_matrix) {
    /* 3x3 identity: best is the diagonal, score 3 */
    std::vector<std::vector<double> > m;
    for (int i = 0; i < 3; ++i) {
        std::vector<double> row(3, 0.0);
        row[i] = 1.0;
        m.push_back(row);
    }
    EnhancedGS egs(m);
    AlignInfo* p = egs.getPAliInfo();
    CHECK_CLOSE(p->sco, 3.0, 1e-9);
    CHECK(p->getAliPairs().size() == 3);
}

TEST(enhancedgs_more_rows_than_cols) {
    /* 3 rows, 2 cols -> at most 2 aligned pairs */
    std::vector<std::vector<double> > m;
    std::vector<double> r0; r0.push_back(0.5); r0.push_back(0.1);
    std::vector<double> r1; r1.push_back(0.1); r1.push_back(0.6);
    std::vector<double> r2; r2.push_back(0.2); r2.push_back(0.2);
    m.push_back(r0); m.push_back(r1); m.push_back(r2);

    EnhancedGS egs(m);
    AlignInfo* p = egs.getPAliInfo();
    CHECK(p->getAliPairs().size() == 2);
    CHECK_CLOSE(p->sco, 1.1, 1e-9);
}
