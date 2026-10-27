#!/usr/bin/env bash
# Black-box: auxiliary output formats and mode flags.
#   -o  superposed query as PDB
#   -m  rotation matrix + transform code
#   -H  load hydrogens
#   -rf 1  flexible mode smoke test (deterministic on the small octane fixture)
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

TMPD="$(mktemp -d)"

echo "== CLI: -o superposition PDB =="
run_ls "${FIX}/self.mol2" "${FIX}/copy.mol2" -rf 0 -o "${TMPD}/sup.pdb" >/dev/null
if [ -s "${TMPD}/sup.pdb" ] && grep -q "^HETATM" "${TMPD}/sup.pdb"; then
    CLI_PASS=$((CLI_PASS+1)); echo "[  OK  ] -o wrote HETATM records"
else CLI_FAIL=$((CLI_FAIL+1)); echo "[ FAIL] -o output missing/short"; fi

echo "== CLI: -m rotation matrix =="
run_ls "${FIX}/self.mol2" "${FIX}/copy.mol2" -rf 0 -m "${TMPD}/mtx.txt" >/dev/null
if [ -s "${TMPD}/mtx.txt" ] && grep -q "u\[i\]\[0\]" "${TMPD}/mtx.txt"; then
    CLI_PASS=$((CLI_PASS+1)); echo "[  OK  ] -m wrote rotation matrix header"
else CLI_FAIL=$((CLI_FAIL+1)); echo "[ FAIL] -m output missing/unexpected"; fi

echo "== CLI: -H loads hydrogens (xtal-lig 59 atoms) =="
# template is the 2nd input; its size is the last column of the first data row.
TS_H=$(run_ls -H 1 "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2" -rf 0 | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/{print $10; exit}')
assert_eq "-H TEMPL_SIZE" "$TS_H" "59"
TS_noH=$(run_ls "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2" -rf 0 | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/{print $10; exit}')
assert_eq "-noH TEMPL_SIZE" "$TS_noH" "29"

echo "== CLI: flexible mode (-rf 1) smoke on octane =="
R1=$(run_ls -rf 1 "${FIX}/octane.mol2" "${FIX}/octane.mol2" | awk 'NR>2 && $3 ~ /^[0-9.]+$/{print $3; exit}')
R0=$(run_ls -rf 0 "${FIX}/octane.mol2" "${FIX}/octane.mol2" | awk 'NR>2 && $3 ~ /^[0-9.]+$/{print $3; exit}')
assert_close "octane self PC-score (rigid)" "$R0" 1.0 1e-3
assert_inrange "octane self PC-score (flexible)" "$R1" 0.0 1.0

rm -rf "${TMPD}"
finish
