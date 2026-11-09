#include "mini_test.h"
#include "FragmentLSalign.h"
#include <vector>
#include <string>
#include <cmath>

/* FragmentLSalign is experimental and not wired into main() (see the file-head
 * warning in FragmentLSalign.h).  T-044 only locks the stale-bond crash on the
 * constructor path; it does not take on the class's known B14 defects. */

/* A ligand whose bondInfoes row survives a setPoints() shrink (same shape as
 * t044MakeStaleBond in test_rotamer.cpp; kept local so the two files do not
 * depend on each other's statics). */
static Ligand t044FragmentStaleLigand(int nAtoms, int aInd, int bInd, int shrinkTo) {
    std::vector<Point3D> pts;
    std::vector<std::string> types;
    for (int i = 0; i < nAtoms; ++i) {
        pts.push_back(Point3D(i * 1.5, 0.0, 0.0));
        types.push_back("C");
    }
    BondInfo b;
    b.aInd = aInd; b.bInd = bInd; b.bt = SINGLE;
    b.bondDis = 1.5; b.bondDisPow2 = 2.25;
    std::vector<BondInfo> bonds; bonds.push_back(b);
    Ligand L("T044F", pts, types, bonds);
    L.setPoints(std::vector<Point3D>(shrinkTo, Point3D(0.0, 0.0, 0.0)));
    return L;
}

/* T-044 sub-step D: FragmentLSalign::shortestRoadBetweenTwoAtoms() built its
 * adjacency matrix from the raw aInd/bInd as well (measured SIGSEGV at
 * FragmentLSalign.h:510 via ctor -> align() -> divideLigandToFragment()).
 * The stale bond must be skipped, and the degenerate 2-atom ligand must still
 * go through the fallback path without crashing.  The score is only required
 * to be finite: FragmentLSalign's own semantics (B14) are out of scope. */
TEST(t044_fragment_ctor_skips_stale_bond) {
    Ligand q = t044FragmentStaleLigand(100, 0, 99, 2);
    Ligand t = t044FragmentStaleLigand(100, 0, 99, 2);
    CHECK(q.size() == 2);
    CHECK(q.getBondInfoes().size() == 1);

    FragmentLSalign f(q, t, 1.0);
    double score = f.getLSscore();
    CHECK(score == score);            /* not NaN */
    CHECK(score >= 0.0);
    CHECK(score < 1e6);
}
