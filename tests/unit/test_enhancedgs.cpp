#include "mini_test.h"
#include "EnhancedGS.h"
#include <vector>
#include <limits>

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

/* F2: an undefined (NaN) cell must never be selected.  The greedy pass keeps
 * its "no score worse than 0" invariant, so the finite cells still win. */
TEST(f2_enhancedgs_ignores_nan_cell) {
    double nanVal = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double> > m;
    std::vector<double> r0; r0.push_back(nanVal); r0.push_back(0.4);
    std::vector<double> r1; r1.push_back(0.3); r1.push_back(nanVal);
    m.push_back(r0); m.push_back(r1);

    EnhancedGS egs(m);
    AlignInfo* p = egs.getPAliInfo();
    CHECK(p->sco == p->sco);
    CHECK(p->sco >= 0.0);
    /* only the two finite cells can be picked: (0,1)=0.4 and (1,0)=0.3 */
    CHECK_CLOSE(p->sco, 0.7, 1e-9);
    std::vector<AlignedPair> pairs = p->getAliPairs();
    CHECK(pairs.size() == 2);
    for (int i = 0; i < (int)pairs.size(); ++i) {
        CHECK(p->getIthAliVec(pairs[i].rowInd) != -1);
    }
}

/* F2: an all-NaN matrix yields an empty, finite-scored alignment. */
TEST(f2_enhancedgs_all_nan_is_empty) {
    double nanVal = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double> > m;
    std::vector<double> r0(2, nanVal);
    std::vector<double> r1(2, nanVal);
    m.push_back(r0); m.push_back(r1);

    EnhancedGS egs(m);
    AlignInfo* p = egs.getPAliInfo();
    CHECK(p->sco == 0.0);
    CHECK(p->getAliPairs().empty());
}
