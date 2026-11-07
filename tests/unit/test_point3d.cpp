#include "mini_test.h"
#include "Point3D.h"
#include "basic_fun.h"
#include <cmath>

TEST(point3d_constructor_and_getters) {
    Point3D p(1.0, 2.0, 3.0);
    CHECK_CLOSE(p.getX(), 1.0, 1e-12);
    CHECK_CLOSE(p.getY(), 2.0, 1e-12);
    CHECK_CLOSE(p.getZ(), 3.0, 1e-12);
    CHECK_CLOSE(p[0], 1.0, 1e-12);
    CHECK_CLOSE(p[1], 2.0, 1e-12);
    CHECK_CLOSE(p[2], 3.0, 1e-12);
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

/* F2: an out-of-range index used to return an uninitialised stack value.  It
 * is now deterministically NaN, and the accessor leaves the object alone. */
TEST(f2_oob_nan) {
    Point3D p(1.0, 2.0, 3.0);
    CHECK(!(p[3] >= 0.0));
    CHECK(!(p[-1] >= 0.0));
    CHECK(!(p[100] >= 0.0));
    CHECK(p[3] != p[3]);
    /* in-range indices are unchanged */
    CHECK(p[0] == 1.0);
    CHECK(p[1] == 2.0);
    CHECK(p[2] == 3.0);

    const Point3D q(4.0, 5.0, 6.0);
    CHECK(q[2] == 6.0);
    CHECK(!(q[3] >= 0.0));
}

/* F1: the three operators are value-semantic; no heap object is handed out. */
TEST(f1_value_semantics) {
    Point3D a(1.5, -2.25, 3.0);
    Point3D b(0.5, 0.75, -1.0);

    Point3D sum = a + b;
    CHECK_CLOSE(sum.getX(), 2.0, 1e-15);
    CHECK_CLOSE(sum.getY(), -1.5, 1e-15);
    CHECK_CLOSE(sum.getZ(), 2.0, 1e-15);

    Point3D diff = a - b;
    CHECK_CLOSE(diff.getX(), 1.0, 1e-15);
    CHECK_CLOSE(diff.getY(), -3.0, 1e-15);
    CHECK_CLOSE(diff.getZ(), 4.0, 1e-15);

    Point3D half = a / 3;
    CHECK_CLOSE(half.getX(), 0.5, 1e-15);
    CHECK_CLOSE(half.getY(), -0.75, 1e-15);
    CHECK_CLOSE(half.getZ(), 1.0, 1e-15);

    /* chained use must compile and stay independent of the operands */
    Point3D chained = (a - b) + b;
    CHECK_CLOSE(chained.getX(), a.getX(), 1e-15);
    CHECK_CLOSE(chained.getY(), a.getY(), 1e-15);
    CHECK_CLOSE(chained.getZ(), a.getZ(), 1e-15);
    CHECK_CLOSE(a.getX(), 1.5, 1e-15);
}

/* F1: a zero-length axis used to divide by zero and yield a NaN matrix.  It
 * must now degrade to the identity rotation (no rotation at all). */
TEST(f1_zero_axis) {
    Point3D zeroAxis(0.0, 0.0, 0.0);
    double** m = generateRotatedCoordinate(zeroAxis, 1.234);

    int h = 0, d = 0;
    bool allFinite = true;
    for (h = 0; h < 3; ++h) {
        for (d = 0; d < 3; ++d) {
            if (!(m[h][d] == m[h][d])) allFinite = false;
            CHECK_CLOSE(m[h][d], (h == d) ? 1.0 : 0.0, 1e-15);
        }
    }
    CHECK(allFinite);
    release2DArr(3, m);

    /* the same matrix for an ordinary unit axis must still be a real rotation */
    Point3D axis(1.0, 0.0, 0.0);
    double** r = generateRotatedCoordinate(axis, 3.14159265358979323846 / 2.0);
    CHECK_CLOSE(r[0][0], 1.0, 1e-12);
    CHECK_CLOSE(r[1][1], 0.0, 1e-12);
    CHECK_CLOSE(r[1][2], -1.0, 1e-12);
    CHECK_CLOSE(r[2][1], 1.0, 1e-12);
    CHECK_CLOSE(r[2][2], 0.0, 1e-12);
    release2DArr(3, r);
}
