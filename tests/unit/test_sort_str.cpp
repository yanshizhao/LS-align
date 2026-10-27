#include "mini_test.h"
#include "Sort.h"
#include "Str.h"
#include <string>
#include <vector>

TEST(sort_descent_index_order) {
    std::vector<double> v;
    v.push_back(0.1); v.push_back(3.0); v.push_back(2.0); v.push_back(5.0);
    int* idx = descentInsertSortIndex(v);
    /* descending by value: 5.0(idx3), 3.0(idx1), 2.0(idx2), 0.1(idx0) */
    CHECK(idx[0] == 3);
    CHECK(idx[1] == 1);
    CHECK(idx[2] == 2);
    CHECK(idx[3] == 0);
    delete[] idx;
}

TEST(sort_ascent_index_order) {
    std::vector<double> v;
    v.push_back(0.1); v.push_back(3.0); v.push_back(2.0); v.push_back(5.0);
    int* idx = ascentInsertSortIndex(v);
    CHECK(idx[0] == 0);
    CHECK(idx[1] == 2);
    CHECK(idx[2] == 1);
    CHECK(idx[3] == 3);
    delete[] idx;
}

TEST(sort_descent_topn) {
    std::vector<double> v;
    v.push_back(0.1); v.push_back(3.0); v.push_back(2.0); v.push_back(5.0);
    int* idx = descentInsertSortIndex(v, 2);   /* only top 2, descending */
    CHECK(idx[0] == 3);
    CHECK(idx[1] == 1);
    delete[] idx;
}

TEST(str_split_single_char) {
    std::vector<std::string> r = split("a--cdd-d", '-');
    /* consecutive delimiters collapse, empty tokens dropped */
    CHECK(r.size() == 3);
    if (r.size() == 3) {
        CHECK(r[0] == "a");
        CHECK(r[1] == "cdd");
        CHECK(r[2] == "d");
    }
}

TEST(str_erase_all_and_upper) {
    /* eraseAll returns a NEW string (does not modify in place). */
    std::string s = eraseAll("C.ar", " .1234567890", 12);
    CHECK(s == "Car");
    std::string t = eraseAll("c.ar", ".", 1);   /* drop the dot -> "car" */
    toUpperString(t);
    CHECK(t == "CAR");
}

TEST(str_trim) {
    CHECK(trim("  hi  ") == "hi");
}
