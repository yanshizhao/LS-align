#ifndef MINI_TEST_H
#define MINI_TEST_H

/*
 * mini_test.h -- self-contained C++98 micro test framework for the LS-align
 * test suite.  No third-party dependency, no package manager, deterministic.
 *
 * Usage:
 *     #include "mini_test.h"
 *     TEST(my_case) { CHECK(1 + 1 == 2); CHECK_CLOSE(a, b, 1e-6); }
 *     XFAIL(known_bug) { CHECK(current_wrong_behavior == expected_right_behavior); }
 *
 * Exit status:
 *     0  -> every TEST passed and every XFAIL still fails (expected)
 *     1  -> some TEST failed, OR some XFAIL unexpectedly passed (bug was fixed)
 *
 * An XFAIL that starts passing (XPASS) turns the run red on purpose: it forces
 * the maintainer to promote it to TEST, so a fixed bug can never silently
 * regress back.
 */

#include <cstdio>
#include <cmath>
#include <vector>

namespace minitest {

struct TestCase {
    const char* name;
    void (*fn)();
    int flags;          /* 1 == xfail */
};

inline std::vector<TestCase>& registry()   { static std::vector<TestCase> r; return r; }
inline int& failures()                     { static int f = 0; return f; }
inline int& checks()                       { static int c = 0; return c; }

struct Registrar {
    Registrar(const char* name, void (*fn)(), int flags) {
        TestCase t;
        t.name  = name;
        t.fn    = fn;
        t.flags = flags;
        registry().push_back(t);
    }
};

inline void check(bool ok, const char* file, int line, const char* expr) {
    checks()++;
    if (!ok) {
        failures()++;
        std::printf("      CHECK(%s) failed  at %s:%d\n", expr, file, line);
    }
}

inline void checkClose(double a, double b, double tol,
                       const char* file, int line, const char* expr) {
    checks()++;
    /* NaN-safe: use the negation of <= so NaN counts as a failure. */
    bool ok = !(std::fabs(a - b) > tol);
    if (!ok) {
        failures()++;
        std::printf("      CHECK_CLOSE(%s) failed  at %s:%d  (%.10g vs %.10g, tol %.4g)\n",
                    expr, file, line, a, b, tol);
    }
}

inline int run() {
    int npass = 0, nfail = 0, nxfail = 0, nxpass = 0;
    for (std::size_t i = 0; i < registry().size(); ++i) {
        TestCase& t = registry()[i];
        int before  = failures();
        bool xf     = (t.flags & 1) != 0;
        std::printf("[ RUN  ] %s%s\n", t.name, xf ? "  (xfail)" : "");
        std::fflush(stdout);
        t.fn();
        bool failed = (failures() > before);
        if (xf) {
            if (failed) { nxfail++; std::printf("[XFAIL] %s\n", t.name); }
            else        { nxpass++; std::printf("[XPASS] %s  <-- UNEXPECTED PASS; bug likely fixed, promote XFAIL->TEST\n", t.name); }
        } else {
            if (failed) { nfail++; std::printf("[ FAIL] %s\n", t.name); }
            else        { npass++; std::printf("[  OK  ] %s\n", t.name); }
        }
        std::fflush(stdout);
    }
    std::printf("\n==== summary: %d passed, %d failed, %d xfail, %d xpass  (checks=%d) ====\n",
                npass, nfail, nxfail, nxpass, checks());
    if (nfail > 0 || nxpass > 0) {
        std::printf("RESULT: FAILURE\n");
        return 1;
    }
    std::printf("RESULT: OK\n");
    return 0;
}

} /* namespace minitest */

#define TEST(name)                                                       \
    static void mt_##name();                                             \
    static ::minitest::Registrar mt_reg_##name(#name, &mt_##name, 0);    \
    static void mt_##name()

#define XFAIL(name)                                                      \
    static void mt_##name();                                             \
    static ::minitest::Registrar mt_reg_##name(#name, &mt_##name, 1);    \
    static void mt_##name()

#define CHECK(c)          ::minitest::check((c), __FILE__, __LINE__, #c)
#define CHECK_CLOSE(a,b,tol) ::minitest::checkClose((a), (b), (tol), __FILE__, __LINE__, #a " ~= " #b)

#endif /* MINI_TEST_H */
