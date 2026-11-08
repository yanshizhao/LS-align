#include "mini_test.h"

/* Every std header and project dependency must come first: Ligands.h is pulled
 * in under "private public", and any std header first included after that point
 * trips the libstdc++ __xfer_bufptrs access-level check. */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cctype>

#include "Point3D.h"
#include "Str.h"
#include "DijkstraSRoad.h"

/* This file must be #included BEFORE test_ligands.cpp, which includes
 * Ligands.h without the macro. */
#define private public
#include "Ligands.h"
#undef private

/* F11: resetCaches() must drop atomBondTypes.  The size-invariant sentinel in
 * getAtomBondTypes() cannot detect a stale cache when the atom count is
 * unchanged, so only the clear itself guarantees freshness. */
TEST(f11_atom_bond_types_cleared_on_setpoints) {
    Ligands ligs("fixtures/self.mol2", false);
    Ligand L = ligs[0];
    CHECK(L.size() == 6);
    CHECK(L.atomBondTypes.empty());          /* nothing built yet */

    L.getAtomBondTypes();                    /* warm: builds 6 entries */
    CHECK((int)L.atomBondTypes.size() == 6);

    std::vector<Point3D> moved;              /* same count, new geometry */
    for (int i = 0; i < L.size(); ++i) {
        Point3D p = L.getIthPoints(i);
        moved.push_back(Point3D(p[0] * 10.0, p[1] * 10.0, p[2] * 10.0));
    }
    L.setPoints(moved);

    CHECK(L.atomBondTypes.empty());          /* RED if resetCaches skips the swap */
}

/* F11: setPoints must free the caches BEFORE swapping in the new points: the
 * matrices are deleted row-by-row using points.size().  Reversed order with a
 * longer new list walks past the end of the old matrix and corrupts the heap
 * (measured: exit 127 under -D_GLIBCXX_ASSERTIONS). */
TEST(f11_resetcaches_ordering_with_growing_points) {
    Ligands ligs("fixtures/self.mol2", false);
    Ligand L = ligs[0];
    CHECK(L.size() == 6);
    CHECK(L.getBondDistMtx() != NULL);       /* build the 6x6 cache */

    std::vector<Point3D> nine;
    for (int i = 0; i < 9; ++i)
        nine.push_back(Point3D((double)i, 0.0, 0.0));
    L.setPoints(nine);                       /* 6 -> 9: grows */

    CHECK(L.size() == 9);
    CHECK_CLOSE(L.getIthPoints(8)[0], 8.0, 0.0);
    CHECK_CLOSE(L.getIthPoints(0)[0], 0.0, 0.0);
    /* deliberately does NOT call getBondDistMtx() here: bondInfoes still holds
     * indexes 0..5 (valid), but keeping the assertion set narrow avoids any
     * dependency on T-041's getter-side guards. */
}
