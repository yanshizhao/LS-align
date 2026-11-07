#include "mini_test.h"
#include "leak_probe.h"
#include "Ligands.h"
#include "ShortestRoad.h"
#include "EnhancedGS.h"
#include "Kabsch.h"
#include "basic_fun.h"
#include "Evaluation.h"
#include "AtomFixInfo.h"
#include "Point3D.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cctype>
#include <climits>
#include <time.h>

/* This file must be #included BEFORE test_rigid.cpp / test_rotamer.cpp /
 * test_xfail.cpp.  RigidLSalign.h, LigRotamerDetector.h and FlexibleLSalign.h
 * are first pulled in here, under "private public", so the cases below can
 * reach the private members and helpers they need.  Every header included
 * above the #define (and the std headers) is unaffected. */
#define private public
#include "RigidLSalign.h"
#include "LigRotamerDetector.h"
#include "FlexibleLSalign.h"
#undef private

/* F14: rotMtx/transVec are bare member arrays.  The constructors used to leave
 * them uninitialised; they now call initIdentityTransform() before align().
 * align() overwrites both for every input reachable today, so the initial
 * write is a second line of defence (see the audit note in the ticket).  This
 * case locks the helper itself and the "constructed object is never garbage"
 * invariant. */

static int isIdentityEntry(int h, int d) {
    return (h == d) ? 1 : 0;
}

static void poisonTransform(RigidLSalign& ali) {
    int h = 0, d = 0;
    for (h = 0; h < 3; ++h) {
        ali.transVec[h] = 42.0;
        for (d = 0; d < 3; ++d) {
            ali.rotMtx[h][d] = 42.0;
        }
    }
}

TEST(f14_ctor_identity) {
    Ligands ligs("fixtures/self.mol2", false);

    RigidLSalign twoArg(ligs[0], ligs[0]);
    RigidLSalign threeArg(ligs[0], ligs[0], 2.0);
    RigidLSalign fiveArg(ligs[0], ligs[0], 2.0, 5, 1);
    RigidLSalign* ctors[3];
    ctors[0] = &twoArg;
    ctors[1] = &threeArg;
    ctors[2] = &fiveArg;

    int k = 0, h = 0, d = 0;
    for (k = 0; k < 3; ++k) {
        /* a constructed object never carries poison/garbage in the transform */
        for (h = 0; h < 3; ++h) {
            CHECK(ctors[k]->transVec[h] == ctors[k]->transVec[h]);
            CHECK(std::fabs(ctors[k]->transVec[h]) < 1e6);
            for (d = 0; d < 3; ++d) {
                CHECK(ctors[k]->rotMtx[h][d] == ctors[k]->rotMtx[h][d]);
                CHECK(std::fabs(ctors[k]->rotMtx[h][d]) < 1e6);
            }
        }

        /* from any state, the helper restores the exact identity transform */
        poisonTransform(*ctors[k]);
        ctors[k]->initIdentityTransform();
        for (h = 0; h < 3; ++h) {
            CHECK(ctors[k]->transVec[h] == 0.0);
            for (d = 0; d < 3; ++d) {
                CHECK(ctors[k]->rotMtx[h][d] == (double)isIdentityEntry(h, d));
            }
        }
    }

    /* every seed degenerates on a zero-atom query: the documented result is the
     * identity transform (the exact fixed value from the F14 behaviour table) */
    Ligands empty("fixtures/empty_ligand.mol2", false);
    CHECK(empty.size() == 1);
    RigidLSalign emptyAli(empty[0], empty[0], 2.0);
    for (h = 0; h < 3; ++h) {
        CHECK(emptyAli.transVec[h] == 0.0);
        for (d = 0; d < 3; ++d) {
            CHECK(emptyAli.rotMtx[h][d] == (double)isIdentityEntry(h, d));
        }
    }
}

/* F4-1: an AlignInfo with no aligned pairs used to be run through Kabsch with
 * empty point sets, which scored sco=57.3 (measured before the fix) and left
 * the caller's u/t untouched.  It must now return a non-NULL empty AlignInfo
 * (sco=0, len=0) and write u=I, t=0. */
TEST(f4_empty_ali_returns_empty_aligninfo) {
    Ligands ligs("fixtures/self.mol2", false);
    RigidLSalign ali(ligs[0], ligs[0], 2.0);

    int h = 0, d = 0;
    AlignInfo zero(0.0, NULL, 0);
    double u[3][3];
    double t[3];
    for (h = 0; h < 3; ++h) {
        t[h] = -777.0;
        for (d = 0; d < 3; ++d) {
            u[h][d] = -777.0;
        }
    }

    AlignInfo* r = ali.nonsequentialUpdateAliInfo(&zero, 5, 3, true, u, t);
    CHECK(r != NULL);
    CHECK(r->sco == 0.0);
    CHECK(r->len == 0);
    CHECK(r->getAliPairs().empty());
    for (h = 0; h < 3; ++h) {
        CHECK(t[h] == 0.0);
        for (d = 0; d < 3; ++d) {
            CHECK(u[h][d] == (double)isIdentityEntry(h, d));
        }
    }
    delete r;

    /* isNeedRotMtxAndTransVec=false: the caller's arrays must stay untouched */
    AlignInfo zero2(0.0, NULL, 0);
    double u2[3][3];
    double t2[3];
    for (h = 0; h < 3; ++h) {
        t2[h] = -777.0;
        for (d = 0; d < 3; ++d) {
            u2[h][d] = -777.0;
        }
    }
    AlignInfo* r2 = ali.nonsequentialUpdateAliInfo(&zero2, 5, 3, false, u2, t2);
    CHECK(r2 != NULL);
    CHECK(r2->sco == 0.0);
    CHECK(r2->len == 0);
    for (h = 0; h < 3; ++h) {
        CHECK(t2[h] == -777.0);
        for (d = 0; d < 3; ++d) {
            CHECK(u2[h][d] == -777.0);
        }
    }
    delete r2;
}

/* F4-2: when the initial rigid solution has no aligned pairs, selecting
 * rotamers from it is meaningless.  Pre-fix this path read past the end of a
 * zero-length sort-index array (measured: SIGSEGV with a 26-rotamer ligand at
 * sn=1) and returned arbitrary rotamers otherwise.  The selection must now be
 * empty (flexible mode degrades to the initial rigid result, i.e. -rf 1
 * behaves like -rf 0 on such input). */
TEST(f4_select_rotamers_empty_init_returns_empty) {
    Ligands q("fixtures/octane.mol2", false);
    Ligands t("fixtures/self.mol2", false);
    Ligands ua("fixtures/unrelatedA.mol2", false);
    Ligands ub("fixtures/unrelatedB.mol2", false);

    FlexibleLSalign f(q[0], t[0], 2.0, false, false, false);
    CHECK(f.queryRotamers.size() == 2);
    CHECK(f.p_initalRLSalign->getFinalAlignedPair().size() > 0);

    /* swap in an initial solution whose alignment is empty */
    delete f.p_initalRLSalign;
    f.p_initalRLSalign = new RigidLSalign(ua[0], ub[0], 0.1);
    CHECK(f.p_initalRLSalign->getFinalAlignedPair().empty());

    /* 2 rotamers > sn, so the initAlign-driven path is taken */
    std::vector<Ligand> sel = f.selectRotamers(1);
    CHECK(sel.empty());
}

/* F4: the empty-alignment early return allocates a new AlignInfo; the caller
 * deletes it through the normal chain.  After construction + one degenerate
 * call + destruction there must be no net live allocation. */
TEST(f4_degenerate_alignment_frees_all) {
    Ligands ligs("fixtures/self.mol2", false);
    { LeakDelta warm; { RigidLSalign w(ligs[0], ligs[0], 2.0); (void)w.getLSscore8QuerySize(); } (void)warm; }

    LeakDelta d;
    {
        RigidLSalign ali(ligs[0], ligs[0], 2.0);
        (void)ali.getLSscore8QuerySize();

        AlignInfo zero(0.0, NULL, 0);
        double u[3][3];
        double t[3];
        AlignInfo* r = ali.nonsequentialUpdateAliInfo(&zero, 5, 3, true, u, t);
        CHECK(r != NULL);
        CHECK(r->sco == 0.0);
        delete r;
    }
    CHECK(d.delta() == 0);
}

/* F2: parseFinalAliInfo used to add a NaN pair score to `sco` *before* the
 * "is this pair real?" test, and NaN compares false, so a garbage pair was
 * also recorded as an aligned pair.  The geometry check must come first:
 * such a pair contributes nothing, is not recorded and is not counted. */
TEST(f2_parse_final_ali_finite) {
    Ligands ligs("fixtures/self.mol2", false);
    RigidLSalign ali(ligs[0], ligs[0], 2.0);

    /* self vs self: all six atoms are aligned at distance 0 */
    int row = 0, col = 0;
    std::vector<AlignedPair> before = ali.getFinalAlignedPair();
    int goodPairs = (int)before.size();
    CHECK(goodPairs > 0);

    for (row = 0; row < 3; ++row) {
        for (col = 0; col < 3; ++col) {
            ali.rotMtx[row][col] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    ali.finalAlignedPair.clear();
    ali.finalAlignedPairDist.clear();
    ali.parseFinalAliInfo();

    /* every rotated query point is NaN, so every distance is NaN */
    CHECK(ali.p_AliInfo->sco == ali.p_AliInfo->sco);
    CHECK(ali.p_AliInfo->sco >= 0.0);
    CHECK(ali.finalAlignedPair.empty());
    CHECK(ali.finalAlignedPairDist.empty());
    CHECK(ali.alisco == ali.alisco);
    CHECK(ali.alisco == 0.0);
}

/* F2: LSscoreFunc must keep the full query x templ matrix even when a pair's
 * distance is undefined; the distance term degrades to 0, never NaN. */
TEST(f2_lsscore_mtx_no_nan) {
    Ligands ligs("fixtures/self.mol2", false);
    Ligands tiny("fixtures/tiny.mol2", false);
    RigidLSalign ali(tiny[0], ligs[0], 2.0);

    double nanVal = std::numeric_limits<double>::quiet_NaN();
    std::vector<Point3D> queryPts;
    queryPts.push_back(Point3D(nanVal, nanVal, nanVal));
    for (int i = 1; i < (int)tiny[0].size(); ++i) {
        queryPts.push_back(tiny[0].getIthPoints(i));
    }

    AlignInfo* p = ali.LSscoreFunc(queryPts, ligs[0].getPoints());
    CHECK(p != NULL);
    CHECK(p->sco == p->sco);
    CHECK(p->sco >= 0.0);
    /* dimensions are dictated by the ligand sizes, not by the guard */
    CHECK(p->len == (int)tiny[0].size());
    /* the distance term degrades to 0, so every row still has usable
     * (mass / bond) scores and the full matrix stays assignable */
    CHECK((int)p->getAliPairs().size() == (int)tiny[0].size());
    delete p;

    /* all-NaN query: the result is still a well-formed, finite-scored info */
    std::vector<Point3D> allNan;
    for (int i = 0; i < (int)tiny[0].size(); ++i) {
        allNan.push_back(Point3D(nanVal, nanVal, nanVal));
    }
    AlignInfo* q = ali.LSscoreFunc(allNan, ligs[0].getPoints());
    CHECK(q != NULL);
    CHECK(q->sco == q->sco);
    CHECK(q->sco >= 0.0);
    CHECK(q->len == (int)tiny[0].size());
    CHECK((int)q->getAliPairs().size() == (int)tiny[0].size());
    delete q;
}

/* F2: the three-point fingerprint helpers must stay finite when the distance
 * matrix contains an undefined entry, and must still find the real extremes
 * among the good entries. */
TEST(f2_three_special_points_no_nan) {
    Ligands q("fixtures/octane.mol2", false);
    Ligands t("fixtures/self.mol2", false);
    FlexibleLSalign f(q[0], t[0], 2.0, false, false, false);

    double nanVal = std::numeric_limits<double>::quiet_NaN();
    std::vector<Point3D> pts;
    pts.push_back(Point3D(0.0, 0.0, 0.0));
    pts.push_back(Point3D(nanVal, nanVal, nanVal));
    pts.push_back(Point3D(3.0, 0.0, 0.0));
    pts.push_back(Point3D(0.0, 4.0, 0.0));
    pts.push_back(Point3D(0.0, 0.0, 12.0));
    pts.push_back(Point3D(0.0, 1.0, 0.0));

    double** mtx = f.calculateDisPow2Mtx(pts);
    double* three = f.calculateThreeSpecicalDistPow2((int)pts.size(), mtx);
    CHECK(three[0] == three[0]);
    CHECK(three[1] == three[1]);
    CHECK(three[2] == three[2]);
    CHECK(three[0] >= 0.0);
    CHECK(three[1] >= 0.0);
    CHECK(three[2] >= 0.0);
    /* farthest finite pair is (0,4,0)-(0,0,12): 16 + 144 = 160 */
    CHECK_CLOSE(three[0], 160.0, 1e-9);

    double sco = f.scoreWithThreeSpecialPts((int)pts.size(), mtx,
        three[0], three[1], three[2]);
    CHECK(sco == sco);
    CHECK(sco >= 0.0);
    CHECK(sco < 1e18);

    release2DArr((int)pts.size(), mtx);
    delete[] three;
}

/* P2: with zero final pairs the "empty vector" sentinel re-ran the whole
 * settlement on every getter call.  The flag must make settlement happen
 * once.  Poisoning alisco/sco gives a deterministic observable: a re-parse
 * overwrites both with its own values. */
TEST(p2_zero_pair_parses_once) {
    Ligands a("fixtures/unrelatedA.mol2", false);
    Ligands b("fixtures/unrelatedB.mol2", false);
    RigidLSalign ali(a[0], b[0], 0.1);
    /* unrelatedA/B with d0=0.1: EGS finds pairs, but none within d0 */
    CHECK(ali.m_isFinalParsed);
    CHECK(ali.finalAlignedPair.empty());
    CHECK(ali.finalAlignedPairDist.empty());
    /* warm up any one-time allocations outside the measured scope */
    (void)ali.getJaccardRatio();

    ali.alisco = -7.0;
    ali.p_AliInfo->sco = -7.0;

    /* every sentinel-guarded entry point must trust the flag: a re-parse
     * would overwrite the poisoned alisco/sco */
    CHECK(ali.getFinalAlignedPair().empty());
    CHECK(ali.getFinalAlignedPairDist().empty());
    CHECK(ali.getJaccardRatio() == 0.0);
    CHECK(ali.getAlignedRMSD() == 0.0);
    CHECK(ali.getAliSco() == -7.0);
    CHECK(ali.toStrVec().size() >= 3);
    CHECK(ali.formatAliInfo().empty());
    CHECK(ali.alisco == -7.0);
    CHECK(ali.p_AliInfo->sco == -7.0);

    LeakDelta d;
    for (int k = 0; k < 1000; ++k) {
        (void)ali.getJaccardRatio();
        (void)ali.getAlignedRMSD();
        (void)ali.getAliSco();
    }
    /* correct behaviour: already settled, so nothing recomputes */
    CHECK(ali.m_isFinalParsed);
    CHECK(ali.alisco == -7.0);
    CHECK(ali.p_AliInfo->sco == -7.0);
    CHECK(d.delta() == 0);
}

/* P2: copying a settled object must copy the flag, otherwise the copy
 * re-parses an already-settled (possibly empty) result on first use. */
TEST(p2_copy_preserves_flag) {
    Ligands a("fixtures/unrelatedA.mol2", false);
    Ligands b("fixtures/unrelatedB.mol2", false);
    RigidLSalign ali(a[0], b[0], 0.1);
    CHECK(ali.m_isFinalParsed);
    CHECK(ali.getFinalAlignedPair().empty());

    ali.alisco = -7.0;
    ali.p_AliInfo->sco = -7.0;

    RigidLSalign copy(ali);
    CHECK(copy.m_isFinalParsed);
    CHECK(copy.finalAlignedPair.empty());
    CHECK(copy.alisco == -7.0);
    CHECK(copy.p_AliInfo->sco == -7.0);   /* own AlignInfo, copied score */

    /* the copy trusts its flag too: no re-parse on first getter */
    (void)copy.getJaccardRatio();
    CHECK(copy.alisco == -7.0);
    CHECK(copy.p_AliInfo->sco == -7.0);
}

/* F9: the copy constructor used to copy rough_sliding_step but not
 * detail_sliding_step, leaving the copy with an uninitialised increment for
 * align()'s detail pass.  align() is not reachable on a copy today, so the
 * copy is the only observable place; it must carry the source's value. */
TEST(f9_detail_step_copied) {
    Ligands ligs("fixtures/self.mol2", false);
    RigidLSalign src(ligs[0], ligs[0], 2.0, 3, 7);
    CHECK(src.rough_sliding_step == 3);
    CHECK(src.detail_sliding_step == 7);

    RigidLSalign copy(src);
    CHECK(copy.rough_sliding_step == 3);
    CHECK(copy.detail_sliding_step == 7);

    /* deep copy stays deep: the two AlignInfo objects are distinct */
    CHECK(copy.p_AliInfo != src.p_AliInfo);
}

/* F3: copying an object whose p_AliInfo is NULL used to dereference the NULL
 * pointer (`new AlignInfo(*NULL)`) and crash.  The copy must be a valid object
 * with a NULL p_AliInfo, and every score accessor on it must return a safe
 * value instead of crashing. */
TEST(f3_copy_null_safe) {
    Ligands a("fixtures/unrelatedA.mol2", false);
    Ligands b("fixtures/unrelatedB.mol2", false);
    RigidLSalign src(a[0], b[0], 0.1);
    CHECK(src.p_AliInfo != NULL);
    CHECK(src.getFinalAlignedPair().empty());

    AlignInfo* saved = src.p_AliInfo;
    src.p_AliInfo = NULL;

    RigidLSalign copy(src);
    CHECK(copy.p_AliInfo == NULL);
    CHECK(copy.getLSscore() == 0.0);
    CHECK(copy.getLSscore8QuerySize() == 0.0);
    CHECK(copy.getLSscore8TemplSize() == 0.0);
    CHECK(copy.getFinalAlignedPair().empty());
    CHECK(copy.getFinalAlignedPairDist().empty());
    CHECK(copy.getAlignedRMSD() == 0.0);
    CHECK(copy.getJaccardRatio() == 0.0);
    CHECK(copy.toStrVec().size() == 7);
    CHECK(copy.formatAliInfo().empty());

    src.p_AliInfo = saved;
}

/* F3: with p_AliInfo NULL the parse guard must make every downstream entry
 * point safe.  The object is unusable, so the guard leaves the parsed flag
 * alone (false): no false "already settled" state, no crash, scores are 0. */
TEST(f3_accessors_null_safe) {
    Ligands ligs("fixtures/self.mol2", false);
    RigidLSalign ali(ligs[0], ligs[0], 2.0);

    AlignInfo* saved = ali.p_AliInfo;
    ali.p_AliInfo = NULL;
    ali.m_isFinalParsed = false;
    ali.finalAlignedPair.clear();
    ali.finalAlignedPairDist.clear();
    ali.alisco = -1.0;

    CHECK(ali.getLSscore() == 0.0);
    CHECK(ali.getLSscore8QuerySize() == 0.0);
    CHECK(ali.getLSscore8TemplSize() == 0.0);
    CHECK(ali.getFinalAlignedPair().empty());
    CHECK(ali.getFinalAlignedPairDist().empty());
    CHECK(ali.getAlignedRMSD() == 0.0);
    CHECK(ali.getJaccardRatio() == 0.0);
    CHECK(ali.getAliSco() == -1.0);
    CHECK(ali.toStrVec().size() == 7);
    CHECK(ali.formatAliInfo().empty());

    /* the guard does not claim settlement on an unusable object */
    CHECK(!ali.m_isFinalParsed);

    ali.p_AliInfo = saved;
    CHECK(ali.getLSscore() > 0.0);
}
