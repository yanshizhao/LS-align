#include "mini_test.h"
#include "Point3D.h"

TEST(point3d_constructor_and_getters) {
    Point3D p(1.0, 2.0, 3.0);
    CHECK_CLOSE(p.getX(), 1.0, 1e-12);
    CHECK_CLOSE(p.getY(), 2.0, 1e-12);
    CHECK_CLOSE(p.getZ(), 3.0, 1e-12);
    CHECK_CLOSE(p[0], 1.0, 1e-12);
    CHECK_CLOSE(p[1], 2.0, 1e-12);
    CHECK_CLOSE(p[2], 3.0, 1e-12);
    CHECK_CLOSE(p.get(0), 1.0, 1e-12);
}

TEST(point3d_distance) {
    Point3D a(0.0, 0.0, 0.0);
    Point3D b(3.0, 4.0, 0.0);          /* classic 3-4-5 */
    CHECK_CLOSE(a.dist(b), 5.0, 1e-12);
    CHECK_CLOSE(a.distPow2(b), 25.0, 1e-12);
    CHECK_CLOSE(a.dist(a), 0.0, 1e-12);
}

TEST(point3d_copy_is_independent) {
    Point3D a(1.0, 1.0, 1.0);
    Point3D b(a);
    CHECK_CLOSE(b.getX(), a.getX(), 1e-12);
}
