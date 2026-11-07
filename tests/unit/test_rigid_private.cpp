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
