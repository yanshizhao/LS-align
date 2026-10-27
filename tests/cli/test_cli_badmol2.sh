#!/usr/bin/env bash
# Black-box: malformed / degenerate mol2 inputs must not crash the binary.
#   E: a ligand containing an isolated heavy atom (B16)
#   F: a ligand with zero atoms (B24)
# Both are expected to FAIL (SIGSEGV, exit 139) before their fixes; they turn
# green once cutoff matrix + empty-matrix handling land.
# A/B (B21/B23) are appended by T-009; C (B24 columns) by T-014.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

echo "== CLI: isolated atom does not crash (B16) =="
assert_alive "E: isolated_atom x self" "${FIX}/isolated_atom.mol2" "${FIX}/self.mol2" -rf 0

echo "== CLI: zero-atom ligand does not crash (B24) =="
assert_alive "F: empty_ligand x self" "${FIX}/empty_ligand.mol2" "${FIX}/self.mol2" -rf 0

finish
