#include "mini_test.h"
#include <cmath>

/*
 * test_oracle.cpp -- paper formulas that live inline in LSalign.cpp::main()
 * and therefore cannot be called directly from a unit test.  We replicate them
 * here as the oracle and assert their numeric properties.
 */

/* d0(Nmin) = 0.55 * (Nmin-9)^(1/3) + 0.15, clamped to [0.1, 4.0];
 * Nmin < 9 uses 0.1.  (LSalign.cpp:331-352) */
static double expected_d0(int Nmin) {
    double d0;
    if (Nmin < 9) {
        d0 = 0.1;
    } else {
        d0 = 0.55 * std::pow(Nmin - 9, 1.0/3.0) + 0.15;
    }
    if (d0 < 0.1) d0 = 0.1;
    if (d0 > 4.0) d0 = 4.0;
    return d0;
}

/* Gumbel survival: P = 1 - exp(-exp(-(x-u)/delta))   (LSalign.cpp:383-384) */
static double pvalue(double score, double u, double delta) {
    return 1.0 - std::exp(-std::exp(-((score - u) / delta)));
}

TEST(oracle_d0_formula_and_clamps) {
    CHECK_CLOSE(expected_d0(8), 0.1, 1e-12);       /* below threshold -> 0.1 */
    CHECK_CLOSE(expected_d0(9), 0.15, 1e-12);      /* (9-9)^(1/3)=0 -> 0.15   */
    CHECK_CLOSE(expected_d0(17),
                0.55 * std::pow(8.0, 1.0/3.0) + 0.15, 1e-12);
    /* monotonic non-decreasing in Nmin */
    CHECK(expected_d0(9) <= expected_d0(50));
    CHECK(expected_d0(50) <= expected_d0(200));
    /* upper clamp */
    CHECK(expected_d0(100000000) <= 4.0 + 1e-12);
}

TEST(oracle_pvalue_range_and_monotonic) {
    double u = 0.2, delta = 0.08;
    double p_low  = pvalue(0.0, u, delta);
    double p_mid  = pvalue(0.5, u, delta);
    double p_high = pvalue(1.0, u, delta);
    CHECK(p_low >= 0.0 && p_low <= 1.0);
    CHECK(p_mid >= 0.0 && p_mid <= 1.0);
    CHECK(p_high >= 0.0 && p_high <= 1.0);
    /* higher score -> smaller p-value (more significant) */
    CHECK(p_low > p_mid);
    CHECK(p_mid > p_high);
}

TEST(oracle_gumbel_params_are_finite) {
    /* u_T/delta_T from LSalign.cpp:359-365 for Nq=Nt=29 */
    double Nq = 29.0, Nt = 29.0;
    double u_T     = -0.563358 +  0.152792*std::log(Nq) +  0.147930*std::log(Nt);
    double delta_T =  0.110170 + -0.007444*std::log(Nq) + -0.008468*std::log(Nt);
    CHECK(delta_T > 0.0);
    CHECK(u_T > -1000.0 && u_T < 1000.0);
}
