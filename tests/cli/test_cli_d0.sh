#!/usr/bin/env bash
# Black-box: default d0(Nmin) differential check.
# The d0 formula lives inline in main() (LSalign.cpp:331-352).  We verify it by
# differential testing:
#   (a) for Nmin<9 (tiny 4-atom fixture) the default d0 clamps to 0.1, so the
#       default run must be byte-identical to an explicit -d0 0.1 run;
#   (b) for a known Nmin>=9 pair, the computed formula value, passed via -d0,
#       must reproduce the default run exactly.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

# extract just the data rows (numeric 3rd column) for stable comparison
data_rows() { run_ls "$@" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/' ; }

echo "== CLI: default d0 == -d0 0.1 when Nmin < 9 =="
A="$(data_rows -rf 0 "${FIX}/tiny.mol2" "${FIX}/tiny.mol2")"
B="$(data_rows -rf 0 -d0 0.1 "${FIX}/tiny.mol2" "${FIX}/tiny.mol2")"
if [ "$A" = "$B" ] && [ -n "$A" ]; then CLI_PASS=$((CLI_PASS+1)); echo "[  OK  ] default d0 == 0.1 for Nmin=4";
else CLI_FAIL=$((CLI_FAIL+1)); echo "[ FAIL] default d0 mismatch for tiny"; echo "  default: $A"; echo "  -d0 0.1: $B"; fi

echo "== CLI: default d0 == formula d0 for a known Nmin >= 9 pair =="
# Derive Nmin from the actual run (ZINC01535869 vs xtal-lig): Nmin = min(QS,TS).
read _pcq _pct _pq _pt _jac _rmsd QS TS < <(row_for ZINC01535869 -rf 0 "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2")
NMIN=$(( QS < TS ? QS : TS ))
D0=$(awk -v n="${NMIN}" 'BEGIN{ printf "%.10f", 0.55*(n-9)^(1.0/3.0)+0.15 }')
echo "   QUERY_SIZE=${QS} TEMPL_SIZE=${TS} -> Nmin=${NMIN}, computed d0=${D0}"
A="$(data_rows -rf 0 -oq 1 "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2")"
B="$(data_rows -rf 0 -d0 "${D0}" -oq 1 "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2")"
if [ "$A" = "$B" ] && [ -n "$A" ]; then CLI_PASS=$((CLI_PASS+1)); echo "[  OK  ] default d0 == formula d0 (Nmin=${NMIN})";
else CLI_FAIL=$((CLI_FAIL+1)); echo "[ FAIL] default d0 != formula d0"; echo "  default : $A"; echo "  formula : $B"; fi

finish
