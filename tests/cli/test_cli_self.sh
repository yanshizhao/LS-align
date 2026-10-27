#!/usr/bin/env bash
# Black-box: self-alignment and rigid-transform invariance must score ~1.0.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

echo "== CLI: self vs self =="
read PC8Q PC8T P5Q P5T JAC RMSD QS TS < <(row_for SELF -rf 0 "${FIX}/self.mol2" "${FIX}/self.mol2")
assert_close "self PC-score8Q" "$PC8Q" 1.0 1e-3
assert_close "self PC-score8T" "$PC8T" 1.0 1e-3
assert_close "self Jaccard"    "$JAC"  1.0 1e-3
assert_close "self RMSD_lb"    "$RMSD" 0.0 1e-3
assert_eq    "self QUERY_SIZE" "$QS"   6

echo "== CLI: self vs rigid-transformed copy =="
read PC8Q PC8T P5Q P5T JAC RMSD QS TS < <(row_for SELF -rf 0 "${FIX}/self.mol2" "${FIX}/copy.mol2")
assert_close "copy PC-score8Q" "$PC8Q" 1.0 1e-3
assert_close "copy Jaccard"    "$JAC"  1.0 1e-3
assert_close "copy RMSD_lb"    "$RMSD" 0.0 1e-3

finish
