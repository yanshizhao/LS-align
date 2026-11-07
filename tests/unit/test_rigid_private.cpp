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
    /* R4 fixed the aligned-bond filter: octane vs self removes all 7 bonds, so
     * its flexible set is now empty (measured 2 -> 0).  unrelatedA aligns only
     * 3 atoms, leaving 2 rotamers, which is what this case needs. */
    Ligands t("fixtures/unrelatedA.mol2", false);
    /* R1 changed the 4th seed, so unrelatedA x unrelatedB at d0=0.1 now keeps
     * one final pair.  tiny x unrelatedB still has none (both pre- and
     * post-R1), which is what this case needs. */
    Ligands ua("fixtures/tiny.mol2", false);
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
    Ligands a("fixtures/tiny.mol2", false);
    Ligands b("fixtures/unrelatedB.mol2", false);
    RigidLSalign ali(a[0], b[0], 0.1);
    /* tiny/unrelatedB with d0=0.1: EGS finds pairs, but none within d0 */
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
    Ligands a("fixtures/tiny.mol2", false);
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
    Ligands a("fixtures/tiny.mol2", false);
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

/* R8: min-max normalisation skipped constant matrices and leaked the raw value
 * (self x self measured 20.0).  A constant matrix now falls back to the
 * theoretical coverage match/min(G_A,G_B); two identical molecules cover the
 * full smaller bond-weight sum, so every search entry must be exactly 1.0. */
TEST(r8_constant_matrix_uses_coverage) {
    Ligands s("fixtures/self.mol2", false);
    Ligands c("fixtures/copy.mol2", false);

    RigidLSalign selfAli(s[0], s[0], 2.0);
    CHECK((int)selfAli.junBondInfoScoMtx.size() == (int)s[0].size());
    int i = 0, j = 0;
    for (i = 0; i < (int)selfAli.junBondInfoScoMtx.size(); ++i) {
        for (j = 0; j < (int)selfAli.junBondInfoScoMtx[i].size(); ++j) {
            CHECK_CLOSE(selfAli.junBondInfoScoMtx[i][j], 1.0, 1e-12);
        }
    }
    CHECK_CLOSE(selfAli.junBondInfoScoMtx[0][0], 1.0, 1e-12);
    CHECK_CLOSE(selfAli.junBondInfoScoMtx[5][5], 1.0, 1e-12);
    /* the report matrix (Jaccard) is untouched: still 1.0 */
    CHECK_CLOSE(selfAli.junBondInfoScoMtx_for_final_sco[0][0], 1.0, 1e-12);

    /* rigid copy is the same constant matrix: 20.0 before the fix, 1.0 after */
    RigidLSalign copyAli(s[0], c[0], 2.0);
    CHECK_CLOSE(copyAli.junBondInfoScoMtx[0][0], 1.0, 1e-12);
    CHECK_CLOSE(copyAli.junBondInfoScoMtx[5][5], 1.0, 1e-12);
}

/* R8: the other constant-matrix case is the opposite one.  An aromatic C6 ring
 * against an amide N6 ring shares only the "ring" signal, so EGS matches every
 * bond pair at 0.5: overlap 2.0, G_A = 23.6, G_B = 22.8.  The fallback must
 * yield 2.0/22.8 = 0.0877 (measured 2.0 before the fix, 0.0877 after), not 1.0. */
TEST(r8_dissimilar_constant_matrix_uses_coverage) {
    Ligands ar("fixtures/aromatic_ring.mol2", false);
    Ligands am("fixtures/amide_ring.mol2", false);
    RigidLSalign ali(ar[0], am[0], 2.0);

    CHECK((int)ali.junBondInfoScoMtx.size() == 6);
    int i = 0, j = 0;
    for (i = 0; i < (int)ali.junBondInfoScoMtx.size(); ++i) {
        for (j = 0; j < (int)ali.junBondInfoScoMtx[i].size(); ++j) {
            CHECK_CLOSE(ali.junBondInfoScoMtx[i][j], 2.0 / 22.8, 1e-12);
        }
    }
    CHECK_CLOSE(ali.junBondInfoScoMtx[0][0], 0.0877192982456140, 1e-12);
    CHECK(ali.junBondInfoScoMtx[0][0] < 0.1);   /* not mistaken for "identical" */

    /* report matrix (Jaccard 2.0/(23.6+22.8-2.0)) is a different quantity and
     * must not be rewritten by the search-matrix fallback */
    CHECK_CLOSE(ali.junBondInfoScoMtx_for_final_sco[0][0], 2.0 / 44.4, 1e-12);
}

/* R8: a matrix with real spread keeps the original min-max path, bit for bit
 * (octane x octane measured: [0][0]=1/3, [0][1]=0, [1][1]=1). */
TEST(r8_varying_matrix_still_minmax) {
    Ligands o("fixtures/octane.mol2", false);
    RigidLSalign ali(o[0], o[0], 2.0);

    CHECK_CLOSE(ali.junBondInfoScoMtx[0][0], 1.0 / 3.0, 1e-12);
    CHECK_CLOSE(ali.junBondInfoScoMtx[0][1], 0.0, 1e-12);
    CHECK_CLOSE(ali.junBondInfoScoMtx[1][1], 1.0, 1e-12);

    int i = 0, j = 0;
    for (i = 0; i < (int)ali.junBondInfoScoMtx.size(); ++i) {
        for (j = 0; j < (int)ali.junBondInfoScoMtx[i].size(); ++j) {
            CHECK(ali.junBondInfoScoMtx[i][j] >= 0.0);
            CHECK(ali.junBondInfoScoMtx[i][j] <= 1.0);
        }
    }
    /* the report matrix is the unchanged Jaccard matrix (measured [0][0]=1) */
    CHECK_CLOSE(ali.junBondInfoScoMtx_for_final_sco[0][0], 1.0, 1e-12);
}

/* R8: min(G_A,G_B) == 0 (neither atom carries a bond) is defined as coverage
 * 1.0 -- "both have no bonds" is consistent, not zero similarity. */
TEST(r8_zero_bond_coverage_is_one) {
    Ligands s("fixtures/single_atom_templ.mol2", false);
    RigidLSalign ali(s[0], s[0], 2.0);

    CHECK((int)ali.junBondInfoScoMtx.size() == 1);
    CHECK_CLOSE(ali.junBondInfoScoMtx[0][0], 1.0, 1e-12);
    /* the report matrix has its own zero guard and stays 0.0 there */
    CHECK_CLOSE(ali.junBondInfoScoMtx_for_final_sco[0][0], 0.0, 1e-12);
}

/* ---- R1: the 4th initial seed block must use junBondInfoScoMtx ---------- *
 * The block used to receive bondInfoScoMtx, the very same matrix as block 3.
 * Both runs are deterministic and the tie test is a strict '<', so the 4th
 * block could never win: it was dead code, and the only ring-aware seed was
 * never used.  The cases below reconstruct align()'s tournament for the last
 * block twice (jun = fixed, bond = pre-fix) and pin the public result. */

struct r1SeedRun {
    AlignInfo* ali;
    double u[3][3];
    double t[3];
};

/* Runs one seed block exactly as align() does: EGS init, then
 * nonsequentialIterAlign(..., rough_sliding_step, 3, ..., false, 10). */
static r1SeedRun r1RunSeed(RigidLSalign& r, const vector<vector<double> >& mtx) {
    r1SeedRun run;
    EnhancedGS egs(mtx);
    AlignInfo* p_init = egs.getPAliInfo();      /* owned by egs */
    p_init->sco = 0.0;

    double* t = new double[3];
    run.ali = r.nonsequentialIterAlign(p_init, r.rough_sliding_step, 3, true, run.u, t, false, 10);
    for (int h = 0; h < 3; ++h) {
        run.t[h] = t[h];
    }
    delete[] t;
    return run;
}

/* Replays align()'s seed tournament + detail refinement with `fourth` as the
 * 4th block, on `scratch` (a RigidLSalign used only as a result container),
 * and returns getLSscore8QuerySize(). */
static double r1ReplayWithFourth(RigidLSalign& scratch, RigidLSalign& src,
                                 const vector<vector<double> >& fourth) {
    const vector<vector<double> >* mtxs[4];
    mtxs[0] = &src.massDiffScoMtx;
    mtxs[1] = &src.vdwRDiffScoMtx;
    mtxs[2] = &src.bondInfoScoMtx;
    mtxs[3] = &fourth;

    double bestSco = -999999.0;
    r1SeedRun best;
    best.ali = NULL;
    for (int b = 0; b < 4; ++b) {
        r1SeedRun run = r1RunSeed(src, *mtxs[b]);
        if (bestSco < run.ali->sco) {
            bestSco = run.ali->sco;
            if (NULL != best.ali) {
                delete best.ali;
                best.ali = NULL;
            }
            best = run;
        } else {
            delete run.ali;
            run.ali = NULL;
        }
    }

    if (NULL != scratch.p_AliInfo) {
        delete scratch.p_AliInfo;
        scratch.p_AliInfo = NULL;
    }
    double u[3][3];
    double t[3];
    int h = 0, d = 0;
    for (h = 0; h < 3; ++h) {
        t[h] = best.t[h];
        for (d = 0; d < 3; ++d) {
            u[h][d] = best.u[h][d];
        }
    }
    scratch.p_AliInfo = src.nonsequentialIterAlign(best.ali, src.detail_sliding_step, 6, true, u, t, true, 1);
    for (h = 0; h < 3; ++h) {
        scratch.transVec[h] = t[h];
        for (d = 0; d < 3; ++d) {
            scratch.rotMtx[h][d] = u[h][d];
        }
    }
    scratch.m_isFinalParsed = false;
    scratch.finalAlignedPair.clear();
    scratch.finalAlignedPairDist.clear();
    scratch.parseFinalAliInfo();

    delete best.ali;
    best.ali = NULL;
    return scratch.getLSscore8QuerySize();
}

/* ZINC03814184 vs the ACE template at d0=2.0: the two matrices are numerically
 * far apart, jun's seed score strictly beats all three other blocks, and the
 * detail-refined result differs from the pre-fix (bond in block 4) result.
 * NOTE (R3): the absolute seed scores below were re-pinned after R3 enabled the
 * multiscale search in the rough phase; the pre-R3 anchors (ZINC03814177) no
 * longer have a strict jun win under the wider search. */
TEST(r1_fourth_seed_uses_jun_matrix) {
    Ligands q("../example/ace_ligands.mol2", false);
    Ligands t("../example/xtal-lig.mol2", false);
    int idx = -1;
    for (int i = 0; i < (int)q.size(); ++i) {
        if (q[i].getLigname() == "ZINC03814184") { idx = i; break; }
    }
    CHECK(idx >= 0);
    if (idx < 0) return;

    RigidLSalign ali(q[idx], t[0], 2.0);
    RigidLSalign scratch(q[idx], t[0], 2.0);

    /* the two matrices really are different, so feeding the same one twice is
     * observable (measured: jun [0][0]=0.0010, bond [0][0]=0.4964) */
    CHECK(ali.junBondInfoScoMtx[0][0] != ali.bondInfoScoMtx[0][0]);
    CHECK(ali.junBondInfoScoMtx[0][0] < 0.01);
    CHECK(ali.bondInfoScoMtx[0][0] > 0.4);

    /* measured seed scores: mass 12.677449, vdW 12.677449, bond 12.507831,
     * jun 14.800251.  jun is the strict winner, so block 4 is selected. */
    r1SeedRun runMass = r1RunSeed(ali, ali.massDiffScoMtx);
    r1SeedRun runVdwr = r1RunSeed(ali, ali.vdwRDiffScoMtx);
    r1SeedRun runBond = r1RunSeed(ali, ali.bondInfoScoMtx);
    r1SeedRun runJun  = r1RunSeed(ali, ali.junBondInfoScoMtx);
    CHECK_CLOSE(runMass.ali->sco, 12.677449, 1e-5);
    CHECK_CLOSE(runVdwr.ali->sco, 12.677449, 1e-5);
    CHECK_CLOSE(runBond.ali->sco, 12.507831, 1e-5);
    CHECK_CLOSE(runJun.ali->sco, 14.800251, 1e-5);
    CHECK(runJun.ali->sco > runMass.ali->sco);
    CHECK(runJun.ali->sco > runVdwr.ali->sco);
    CHECK(runJun.ali->sco > runBond.ali->sco);
    delete runMass.ali; delete runVdwr.ali; delete runBond.ali; delete runJun.ali;

    /* behaviour: the shipped align() result must equal the tournament whose
     * 4th block is junBondInfoScoMtx, and must differ from the pre-fix
     * tournament whose 4th block was bondInfoScoMtx.  If align() regresses to
     * bondInfoScoMtx in block 4 the first CHECK fails (measured 0.547940 vs
     * 0.649551); if the 4th block is dead code the second CHECK fails. */
    double withJun  = r1ReplayWithFourth(scratch, ali, ali.junBondInfoScoMtx);
    double withBond = r1ReplayWithFourth(scratch, ali, ali.bondInfoScoMtx);

    CHECK_CLOSE(ali.getLSscore8QuerySize(), 0.649550537410, 1e-9);
    CHECK_CLOSE(withJun, 0.649550537410, 1e-9);
    CHECK_CLOSE(withBond, 0.547939632780, 1e-9);
    CHECK_CLOSE(ali.getLSscore8QuerySize(), withJun, 1e-12);
    CHECK(std::fabs(ali.getLSscore8QuerySize() - withBond) > 1e-3);
}

/* ZINC03814200 with the default d0 the CLI would pick: Nmin = 19, so
 * d0 = 0.55*(19-9)^(1/3) + 0.15 = 1.334939...  With the fixed code the
 * final score is 0.6700748; the pre-fix block-4 matrix lands on 0.6237836.
 * This pins the default-d0 path the CLI actually uses.
 * NOTE (R3): re-anchored from ZINC03814164 (d0=1.25 exactly) because under
 * R3's multiscale rough search its 4th block no longer changes the result
 * (jun and bond replays both land on the shipped score). */
TEST(r1_fourth_seed_differs_from_bond_on_ace) {
    Ligands q("../example/ace_ligands.mol2", false);
    Ligands t("../example/xtal-lig.mol2", false);
    int idx = -1;
    for (int i = 0; i < (int)q.size(); ++i) {
        if (q[i].getLigname() == "ZINC03814200") { idx = i; break; }
    }
    CHECK(idx >= 0);
    if (idx < 0) return;

    int nMin = (int)q[idx].size() < (int)t[0].size() ? (int)q[idx].size() : (int)t[0].size();
    double d0 = 0.55 * std::pow(nMin - 9.0, 1.0 / 3.0) + 0.15;
    CHECK_CLOSE(d0, 1.334939079518, 1e-9);

    RigidLSalign ali(q[idx], t[0], d0);
    RigidLSalign scratch(q[idx], t[0], d0);

    double withJun  = r1ReplayWithFourth(scratch, ali, ali.junBondInfoScoMtx);
    double withBond = r1ReplayWithFourth(scratch, ali, ali.bondInfoScoMtx);
    CHECK_CLOSE(ali.getLSscore8QuerySize(), withJun, 1e-12);
    /* measured: 0.6700748 (jun, shipped) vs 0.6237836 (pre-fix) */
    CHECK_CLOSE(ali.getLSscore8QuerySize(), 0.670074793318, 1e-9);
    CHECK_CLOSE(withBond, 0.623783648256, 1e-9);
    CHECK(std::fabs(ali.getLSscore8QuerySize() - withBond) > 1e-3);
}

/* ---- R3: multiscale sliding window ------------------------------------- *
 * The window list was guarded by `if (L_min <= 4)`, and L_min = min(n_ali,4)
 * is always <= 4, so the multiscale `else` branch was dead and every call ran
 * a single window of length L_min.  The guard is now `if (n_ali <= 4)`, so an
 * alignment longer than 4 atoms searches {n_ali, n_ali/2, ..., 4} (capped by
 * max_win_type_num), while n_ali <= 4 is unchanged.
 *
 * L_init_vec is a local of nonsequentialUpdateAliInfo, so the cases observe
 * the real function instead of re-implementing its window loop:
 *   - max_win_type_num = 0 is the pre-R3 list exactly ({L_min});
 *   - more windows can only raise the best score (strict `>` selection);
 *   - each (window, slide) pair allocates a fixed number of blocks, so the
 *     test TU's allocation counter (leak_probe.h) distinguishes one window
 *     from many without depending on the loop's internals.
 * Measured on the 40-atom synthetic pair below, step=1:
 *   mw0=18213  mw1=18714  mw2=29173  mw3=44549  mw4=62334  mw5=62334  mw6=62334
 * The three plateaus are the lists {4}, {40,4}, {40,20,4}, {40,20,10,4} and
 * {40,20,10,5,4} -> widening max_win_type_num past 4 cannot add another
 * entry, which is exactly "5 windows for n_ali=40". */

static AlignInfo* r3Ali40(RigidLSalign& ali, int nPairs) {
    int n = (int)ali.query.size();
    int* alivec = new int[n];
    int i = 0;
    for (i = 0; i < n; ++i) alivec[i] = -1;
    for (i = 0; i < nPairs; ++i) alivec[i] = i;
    AlignInfo* ai = new AlignInfo(0.0, alivec, n);
    delete[] alivec;
    return ai;
}

static long r3AllocsForWindows(RigidLSalign& ali, int nPairs, int maxWin) {
    AlignInfo* ai = r3Ali40(ali, nPairs);
    double u[3][3];
    double t[3];
    long before = totalAllocs();
    AlignInfo* r = ali.nonsequentialUpdateAliInfo(ai, 1, maxWin, false, u, t);
    long after = totalAllocs();
    delete r;
    delete ai;
    return after - before;
}

TEST(r3_window_list_grows_with_n_ali) {
    Ligands q("fixtures/r3_query40.mol2", false);
    Ligands t("fixtures/r3_templ40.mol2", false);
    CHECK(q[0].size() == 40);
    CHECK(t[0].size() == 40);
    RigidLSalign ali(q[0], t[0], 2.0);

    /* n_ali = 40: the single-window list (max_win_type_num 0, = pre-R3) is
     * strictly cheaper than the multiscale list, and the list saturates at
     * max_win_type_num = 4 -> {40,20,10,5,4}, five entries. */
    long single = r3AllocsForWindows(ali, 40, 0);
    long mw1 = r3AllocsForWindows(ali, 40, 1);
    long mw2 = r3AllocsForWindows(ali, 40, 2);
    long mw3 = r3AllocsForWindows(ali, 40, 3);
    long mw4 = r3AllocsForWindows(ali, 40, 4);
    long mw6 = r3AllocsForWindows(ali, 40, 6);
    CHECK(single == 18213);
    CHECK(mw1 == 18714);
    CHECK(mw2 == 29173);
    CHECK(mw3 == 44549);
    CHECK(mw4 == 62334);
    CHECK(mw6 == 62334);          /* cap reached at 4, not at 6 */
    CHECK(mw1 > single);          /* {40,4} adds work */
    CHECK(mw2 > mw1);             /* {40,20,4} */
    CHECK(mw3 > mw2);             /* {40,20,10,4} */
    CHECK(mw4 > mw3);             /* {40,20,10,5,4} */
    CHECK(mw4 == mw6);

    /* n_ali <= 4: always exactly one window ({n_ali}); max_win_type_num has
     * no effect at all, matching the pre-R3 behaviour. */
    int n = 0;
    int mw = 0;
    for (n = 1; n <= 4; ++n) {
        long base = r3AllocsForWindows(ali, n, 0);
        for (mw = 1; mw <= 6; ++mw) {
            CHECK(r3AllocsForWindows(ali, n, mw) == base);
        }
    }
    CHECK(r3AllocsForWindows(ali, 1, 6) == 491);
    CHECK(r3AllocsForWindows(ali, 4, 6) == 497);

    /* n_ali = 9: {9,4} saturates at max_win_type_num = 1, so the single
     * window's count differs from the two-entry list while mw1..mw6 agree.
     * Measured 2959 -> 3456. */
    CHECK(r3AllocsForWindows(ali, 9, 0) == 2959);
    CHECK(r3AllocsForWindows(ali, 9, 1) == 3456);
    CHECK(r3AllocsForWindows(ali, 9, 6) == 3456);

    /* n_ali = 13: {13,6,4}, i.e. exactly two extra windows.  Measured
     * 4927 -> 5424 -> 9377, saturated from max_win_type_num = 2. */
    CHECK(r3AllocsForWindows(ali, 13, 0) == 4927);
    CHECK(r3AllocsForWindows(ali, 13, 1) == 5424);
    CHECK(r3AllocsForWindows(ali, 13, 2) == 9377);
    CHECK(r3AllocsForWindows(ali, 13, 6) == 9377);
}

/* The window search can only take the best score (it uses a strict '>'), so
 * widening the window list must never lower the winner.  Numbers measured on
 * ZINC01535869 x the ACE template at the CLI's default d0. */
TEST(r3_more_windows_not_worse) {
    Ligands q("../example/ace_ligands.mol2", false);
    Ligands t("../example/xtal-lig.mol2", false);
    int idx = -1;
    for (int i = 0; i < (int)q.size(); ++i) {
        if (q[i].getLigname() == "ZINC01535869") { idx = i; break; }
    }
    CHECK(idx >= 0);
    if (idx < 0) return;

    int nMin = (int)q[idx].size() < (int)t[0].size() ? (int)q[idx].size() : (int)t[0].size();
    double d0 = 0.55 * std::pow(nMin - 9.0, 1.0 / 3.0) + 0.15;
    RigidLSalign ali(q[idx], t[0], d0);

    /* a long, valid alignment: the shipped 13-pair result of this pair */
    std::vector<AlignedPair> pairs = ali.getFinalAlignedPair();
    CHECK((int)pairs.size() == 13);
    int n = (int)ali.query.size();
    int* alivec = new int[n];
    int i = 0;
    for (i = 0; i < n; ++i) alivec[i] = -1;
    for (i = 0; i < (int)pairs.size(); ++i) alivec[pairs[i].rowInd] = pairs[i].colInd;
    AlignInfo baseInfo(0.0, alivec, n);
    delete[] alivec;

    AlignInfo single(baseInfo);
    double u[3][3];
    double t3[3];
    AlignInfo* oldR = ali.nonsequentialUpdateAliInfo(&single, 1, 0, false, u, t3);
    CHECK_CLOSE(oldR->sco, 10.207368797088, 1e-9);   /* pre-R3 single window {4} */
    delete oldR;

    AlignInfo multi(baseInfo);
    AlignInfo* newR = ali.nonsequentialUpdateAliInfo(&multi, 1, 6, false, u, t3);
    CHECK_CLOSE(newR->sco, 10.436078061110, 1e-9);   /* {13,6,4} wins */
    CHECK(newR->sco > 10.207368797088);              /* strictly better, never worse */
    delete newR;
}
