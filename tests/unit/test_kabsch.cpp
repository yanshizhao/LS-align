#include "mini_test.h"
#include "Kabsch.h"
#include "Point3D.h"
#include <cmath>
#include <vector>

/* Build a set of non-degenerate points (a 3D "L" shape) so Kabsch is
 * well-conditioned (the two principal axes are distinct). */
static std::vector<Point3D> makeL() {
    std::vector<Point3D> v;
    v.push_back(Point3D(0.0, 0.0, 0.0));
    v.push_back(Point3D(2.0, 0.0, 0.0));
    v.push_back(Point3D(0.0, 3.0, 0.0));
    v.push_back(Point3D(0.0, 0.0, 1.5));
    v.push_back(Point3D(2.0, 3.0, 1.5));
    return v;
}

/* Apply the same transform convention used everywhere in the code:
 *   X = t + u * x
 */
static std::vector<Point3D> applyT(std::vector<Point3D> x,
                                   const double t[3], const double u[3][3]) {
    std::vector<Point3D> out;
    for (std::size_t i = 0; i < x.size(); ++i) {
        double X = t[0] + u[0][0]*x[i][0] + u[0][1]*x[i][1] + u[0][2]*x[i][2];
        double Y = t[1] + u[1][0]*x[i][0] + u[1][1]*x[i][1] + u[1][2]*x[i][2];
        double Z = t[2] + u[2][0]*x[i][0] + u[2][1]*x[i][1] + u[2][2]*x[i][2];
        out.push_back(Point3D(X, Y, Z));
    }
    return out;
}

TEST(kabsch_recovers_known_rotation_about_z) {
    double th = 37.0 * 3.14159265358979323846 / 180.0;
    double R[3][3];
    R[0][0] =  std::cos(th); R[0][1] = -std::sin(th); R[0][2] = 0.0;
    R[1][0] =  std::sin(th); R[1][1] =  std::cos(th); R[1][2] = 0.0;
    R[2][0] =  0.0;          R[2][1] =  0.0;          R[2][2] = 1.0;
    double T[3]; T[0] = 10.0; T[1] = -5.0; T[2] = 3.0;

    std::vector<Point3D> x = makeL();
    std::vector<Point3D> y = applyT(x, T, R);

    double t[3], u[3][3];
    Kabsch(x, y, t, u);

    /* u should reproduce R, t should reproduce T */
    for (int i = 0; i < 3; ++i) {
        CHECK_CLOSE(t[i], T[i], 1e-6);
        for (int j = 0; j < 3; ++j) {
            CHECK_CLOSE(u[i][j], R[i][j], 1e-6);
        }
    }
}

TEST(kabsch_identity_when_sets_equal) {
    std::vector<Point3D> x = makeL();
    double t[3], u[3][3];
    Kabsch(x, x, t, u);
    CHECK_CLOSE(t[0], 0.0, 1e-6);
    CHECK_CLOSE(t[1], 0.0, 1e-6);
    CHECK_CLOSE(t[2], 0.0, 1e-6);
    CHECK_CLOSE(u[0][0], 1.0, 1e-6);
    CHECK_CLOSE(u[1][1], 1.0, 1e-6);
    CHECK_CLOSE(u[2][2], 1.0, 1e-6);
}

TEST(kabsch_overload_with_index_ranges) {
    double th = 25.0 * 3.14159265358979323846 / 180.0;
    double R[3][3];
    R[0][0] = std::cos(th); R[0][1] = -std::sin(th); R[0][2] = 0.0;
    R[1][0] = std::sin(th); R[1][1] =  std::cos(th); R[1][2] = 0.0;
    R[2][0] = 0.0;          R[2][1] =  0.0;          R[2][2] = 1.0;
    double T[3]; T[0] = 1.0; T[1] = 2.0; T[2] = -1.0;

    std::vector<Point3D> x = makeL();
    std::vector<Point3D> y = applyT(x, T, R);

    double t[3], u[3][3];
    Kabsch(x, 0, (int)x.size(), y, 0, (int)y.size(), t, u);
    for (int i = 0; i < 3; ++i) {
        CHECK_CLOSE(t[i], T[i], 1e-6);
        for (int j = 0; j < 3; ++j) CHECK_CLOSE(u[i][j], R[i][j], 1e-6);
    }
}
