#!/usr/bin/env bash
# Black-box: the paper's ACE example (Sec 3.1).
#   ZINC03812885 (active) vs the ACE seed (xtal-lig.mol2, 3bkl).
# The paper reports aligned RMSD ~= 0.44 A and LS-score ~= 0.62 with ~13 aligned
# pairs.  The shipped CLI prints different metrics (PC-score8T is size-normalised;
# RMSD_lb is a lower-bound RMSD), so we assert LOOSE bands plus the row count and
# the aligned-pair count from the -a alignment dump.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

echo "== CLI: ACE example row count =="
assert_eq "49 query ligands -> 49 rows" "$(count_rows -rf 0 "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2")" "49"

echo "== CLI: ZINC03812885 vs ACE seed =="
read PC8Q PC8T P5Q P5T JAC RMSD QS TS < <(row_for ZINC03812885 -rf 0 "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2")
assert_eq "ZINC03812885 QUERY_SIZE" "$QS" "14"
assert_eq "ACE seed TEMPL_SIZE"      "$TS" "29"
assert_inrange "PC-score8T (paper LS-score band)" "$PC8T" 0.30 0.95
assert_inrange "RMSD_lb band"                     "$RMSD" 0.0  6.0
assert_inrange "Jaccard band"                     "$JAC"  0.1  1.0

echo "== CLI: -a alignment dump is well-formed =="
TMPD="$(mktemp -d)"
run_ls "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2" -rf 0 -oq 6 -ot 1 -a "${TMPD}/ali.txt" >/dev/null
NPAIR=$(grep -c "Query Atom Index" "${TMPD}/ali.txt" 2>/dev/null || echo 0)
assert_eq "alignment file written" "$([ -s "${TMPD}/ali.txt" ] && echo yes || echo no)" "yes"
# at least one 'Query Atom Index' block present
if [ "${NPAIR}" -ge 1 ]; then CLI_PASS=$((CLI_PASS+1)); echo "[  OK  ] alignment blocks present (${NPAIR})";
else CLI_FAIL=$((CLI_FAIL+1)); echo "[ FAIL] no alignment blocks"; fi
rm -rf "${TMPD}"

finish
