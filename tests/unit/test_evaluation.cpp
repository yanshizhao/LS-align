#include "mini_test.h"
#include "Evaluation.h"
#include "Point3D.h"
#include <vector>
#include <limits>

TEST(rmsd_lb_identical_sets_is_zero) {
    std::vector<Point3D> a;
    a.push_back(Point3D(0.0, 0.0, 0.0));
    a.push_back(Point3D(10.0, 0.0, 0.0));
    CHECK_CLOSE(RMSD_lb(a, a), 0.0, 1e-9);
}

TEST(rmsd_lb_uniform_translation) {
    std::vector<Point3D> a, b;
    a.push_back(Point3D(0.0, 0.0, 0.0));
    a.push_back(Point3D(10.0, 0.0, 0.0));
    b.push_back(Point3D(1.0, 0.0, 0.0));
    b.push_back(Point3D(11.0, 0.0, 0.0));
    /* every a-atom's nearest b-atom is exactly 1 away, and vice versa */
    CHECK_CLOSE(RMSD_lb(a, b), 1.0, 1e-9);
}

TEST(rmsd_plain_when_sizes_match) {
    std::vector<Point3D> a, b;
    a.push_back(Point3D(0.0, 0.0, 0.0));
    a.push_back(Point3D(3.0, 4.0, 0.0));
    b.push_back(Point3D(0.0, 0.0, 0.0));
    b.push_back(Point3D(0.0, 0.0, 0.0));
    /* sqrt((0 + 25)/2) */
    CHECK_CLOSE(RMSD(a, b), std::sqrt(12.5), 1e-9);
}

/* F2: a NaN distance must not poison the totals or the nearest-neighbour
 * minima; all four helpers return a finite value for the good points. */
static std::vector<Point3D> withNaNPt() {
    std::vector<Point3D> v;
    double nanVal = std::numeric_limits<double>::quiet_NaN();
    v.push_back(Point3D(0.0, 0.0, 0.0));
    v.push_back(Point3D(nanVal, nanVal, nanVal));
    v.push_back(Point3D(2.0, 0.0, 0.0));
    return v;
}

TEST(f2_rmsd_skips_nan_distance) {
    std::vector<Point3D> a = withNaNPt();
    std::vector<Point3D> b = withNaNPt();
    /* finite points contribute 0; the NaN pair is ignored, not summed */
    CHECK(RMSD(a, b) == 0.0);
    CHECK(RMSD_ub(a, b) == 0.0);
}

TEST(f2_rmsd_lb_stays_finite_on_nan) {
    double nanVal = std::numeric_limits<double>::quiet_NaN();
    std::vector<Point3D> a;
    a.push_back(Point3D(0.0, 0.0, 0.0));
    a.push_back(Point3D(nanVal, nanVal, nanVal));

    std::vector<Point3D> b;
    b.push_back(Point3D(1.0, 0.0, 0.0));
    b.push_back(Point3D(3.0, 0.0, 0.0));

    /* the NaN row/column can never win the minimum: it keeps the 1e17 cap */
    double lb = RMSD_lb(a, b);
    double lbMin = RMSD_lb_min(a, b);
    CHECK(lb == lb);
    CHECK(lbMin == lbMin);
    CHECK(lb >= 0.0);
    CHECK(lbMin >= 0.0);
    CHECK(lb <= 1e17);
    CHECK(lbMin <= 1e17);
}
