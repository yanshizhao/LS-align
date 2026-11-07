#!/usr/bin/env bash
# Black-box: strict argument parsing, exit codes, stdout purity (B27/B28).
#
#   * valid runs exit 0; -h / -v exit 0
#   * invalid option values exit 1 with an ERROR on stderr and NO stdout table
#   * -rf 2 names Fragment mode explicitly
#   * library diagnostics live on stderr, never on stdout
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

SCRATCH="$(mktemp -d)"
trap 'rm -rf "${SCRATCH}"' EXIT

stdout_rows() {   # $1 = file with captured stdout
    awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}' "$1"
}

# reject_case LABEL ARGS...  -- expect exit 1, an ERROR on stderr, no stdout rows.
reject_case() {
    local label="$1"; shift
    local out_f="${SCRATCH}/o" err_f="${SCRATCH}/e" rc err
    "${BIN}" "$@" >"${out_f}" 2>"${err_f}"
    rc=$?
    err="$(cat "${err_f}")"
    assert_eq "${label}: exit 1" "${rc}" "1"
    if printf '%s' "${err}" | grep -q "ERROR"; then
        ok "${label}: stderr reports ERROR"
    else
        bad "${label}: stderr missing ERROR (got: ${err})"
    fi
    assert_eq "${label}: no stdout rows" "$(stdout_rows "${out_f}")" "0"
}

QFIX="${FIX}/self.mol2"
TFIX="${FIX}/copy.mol2"

echo "== CLI: exit codes =="
"${BIN}" "${QFIX}" "${TFIX}" -rf 0 >"${SCRATCH}/o" 2>/dev/null
assert_eq "valid -rf 0 run exits 0" "$?" "0"
assert_close "valid run still scores 1.0" \
    "$(awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {print $3; exit}' "${SCRATCH}/o")" 1.0 1e-3
"${BIN}" -h >/dev/null 2>&1
assert_eq "-h exits 0" "$?" "0"
"${BIN}" -v >/dev/null 2>&1
assert_eq "-v exits 0" "$?" "0"

echo "== CLI: invalid -rf =="
reject_case "-rf 3"     "${QFIX}" "${TFIX}" -rf 3
reject_case "-rf -1"    "${QFIX}" "${TFIX}" -rf -1
reject_case "-rf abc"   "${QFIX}" "${TFIX}" -rf abc
reject_case "-rf (no value)" "${QFIX}" "${TFIX}" -rf

echo "== CLI: -rf 2 names the unsupported Fragment mode =="
run_ls_capture rf2_out rf2_err "${QFIX}" "${TFIX}" -rf 2
assert_eq "-rf 2 exit 1" "${RC}" "1"
if printf '%s' "${rf2_err}" | grep -q "Fragment mode"; then
    ok "-rf 2 stderr mentions Fragment mode"
else
    bad "-rf 2 stderr missing 'Fragment mode' (got: ${rf2_err})"
fi
assert_eq "-rf 2: no stdout rows" "$(printf '%s' "${rf2_out}" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}')" "0"

echo "== CLI: invalid 0/1 flags and -d0 =="
reject_case "-md 0.5"  "${QFIX}" "${TFIX}" -md 0.5
reject_case "-md 2"    "${QFIX}" "${TFIX}" -md 2
reject_case "-acc 0.5" "${QFIX}" "${TFIX}" -acc 0.5
reject_case "-H 0.5"   "${QFIX}" "${TFIX}" -H 0.5
reject_case "-d0 abc"  "${QFIX}" "${TFIX}" -d0 abc
reject_case "-d0 0.05" "${QFIX}" "${TFIX}" -d0 0.05
reject_case "-d0 9"    "${QFIX}" "${TFIX}" -d0 9

echo "== CLI: invalid / out-of-range -oq and -ot =="
reject_case "-oq 0"    "${QFIX}" "${TFIX}" -oq 0
reject_case "-oq abc"  "${QFIX}" "${TFIX}" -oq abc
reject_case "-oq 999"  "${QFIX}" "${TFIX}" -oq 999
reject_case "-ot 999"  "${QFIX}" "${TFIX}" -ot 999
# -oq 999 with a real 49-ligand query file must name the count in the error.
run_ls_capture oq_out oq_err -oq 999 "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2" -rf 0
assert_eq "-oq 999 (49 ligands) exit 1" "${RC}" "1"
if printf '%s' "${oq_err}" | grep -q "out of range (query file contains 49 ligand"; then
    ok "-oq 999 stderr states the ligand count"
else
    bad "-oq 999 stderr missing count (got: ${oq_err})"
fi
assert_eq "-oq 999: no stdout rows" "$(printf '%s' "${oq_out}" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}')" "0"

echo "== CLI: in-range -oq/-ot still works =="
"${BIN}" "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2" -oq 6 -ot 1 -rf 0 >"${SCRATCH}/o" 2>/dev/null
assert_eq "-oq 6 -ot 1 exits 0" "$?" "0"
assert_eq "-oq 6 -ot 1 emits one row" "$(stdout_rows "${SCRATCH}/o")" "1"

echo "== CLI: diagnostics stay off stdout (B27) =="
run_ls_capture diag_out diag_err "${FIX}/isolated_atom.mol2" "${FIX}/isolated_atom.mol2" -rf 0
assert_eq "isolated x isolated exits 0" "${RC}" "0"
if printf '%s' "${diag_out}" | grep -q "has a problem"; then
    bad "stdout contains 'has a problem'"
else
    ok "stdout is free of 'has a problem'"
fi
if printf '%s' "${diag_out}" | grep -q "scoMtx\["; then
    bad "stdout contains a scoMtx diagnostic"
else
    ok "stdout is free of scoMtx diagnostics"
fi
# Count is measured, not contractual: it changes when the initial seeds change
# (R1: the junBond 4th seed turns 8 into 4 here).  The contract is only that
# these diagnostics land on stderr, never on stdout.
NDIAG="$(printf '%s\n' "${diag_err}" | grep -c "has a problem")"
assert_eq "stderr carries 4 'has a problem' lines" "${NDIAG}" "4"
# The table must still be well-formed: blank, 2 header lines, 1 data row, blank, 3 timing lines.
BODY_ROWS="$(printf '%s\n' "${diag_out}" | awk 'NF {n++} END{print n+0}')"
assert_eq "stdout has 6 non-blank lines (2 header + 1 row + 3 timing)" "${BODY_ROWS}" "6"
FIRST_AFTER_HEADER="$(printf '%s\n' "${diag_out}" | awk '/\(Without "H"\)/ {getline; print; exit}')"
if printf '%s' "${FIRST_AFTER_HEADER}" | grep -q "ISOLATED"; then
    ok "header is immediately followed by the data row"
else
    bad "data row is not adjacent to the header (got: ${FIRST_AFTER_HEADER})"
fi

finish
