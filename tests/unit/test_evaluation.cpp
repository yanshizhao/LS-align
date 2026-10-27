#include "mini_test.h"
#include "Evaluation.h"
#include "Point3D.h"
#include <vector>

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
