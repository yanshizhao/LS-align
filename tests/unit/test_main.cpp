/*
 * test_main.cpp -- single TU for the white-box unit test suite.
 *
 * It includes the test source files directly (each registers TEST/XFAIL cases
 * into the mini_test registry at static-init time), then runs them.  It does
 * NOT include LSalign.cpp, because that file owns main(); the headers are all
 * header-only, so nothing else is needed.
 *
 * Run from the tests/ directory (fixtures/ paths are relative to it).
 */

#include "mini_test.h"
#include "leak_probe.h"     /* replaces global new/delete in THIS TU only */

#include "test_point3d.cpp"
#include "test_sort_str.cpp"
#include "test_kabsch.cpp"
#include "test_evaluation.cpp"
#include "test_enhancedgs.cpp"
/* Must precede test_ligands: it pulls Ligands.h in under "private public". */
#include "test_ligands_private.cpp"
#include "test_ligands.cpp"
/* Must precede test_rigid/test_rotamer/test_xfail: it pulls RigidLSalign.h,
 * LigRotamerDetector.h and FlexibleLSalign.h in under "private public". */
#include "test_rigid_private.cpp"
#include "test_rigid.cpp"
#include "test_rotamer.cpp"
/* After test_rotamer.cpp: FragmentLSalign.h is first included here, outside the
 * "private public" window above (its constructor is public). */
#include "test_fragment.cpp"
#include "test_oracle.cpp"
#include "test_xfail.cpp"
#include "test_delete_mismatch.cpp"

int main() {
    return ::minitest::run();
}
