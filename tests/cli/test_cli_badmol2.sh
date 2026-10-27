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

# ---------------------------------------------------------------------------
# A/B: malformed mol2 must be rejected with a diagnostic, not crash.
# The binary reports to stderr, so these use run_ls_capture rather than run_ls.
# Line numbers are derived from the fixture at test time (never hard-coded).
# ---------------------------------------------------------------------------
data_rows() {   # $1 = stdout text
    printf '%s\n' "$1" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}'
}

echo "== CLI: truncated ATOM record is rejected (B21) =="
short_line="$(awk '
    /^@<TRIPOS>ATOM/ { inatom=1; next }
    /^@<TRIPOS>/     { inatom=0 }
    inatom && length($0) > 0 && length($0) < 48 { print NR; exit }
' "${FIX}/bad_short_atom.mol2")"
run_ls_capture a_out a_err "${FIX}/bad_short_atom.mol2" "${FIX}/self.mol2" -rf 0
if printf '%s' "${a_err}" | grep -q "ATOM record too short"; then
    ok "A: stderr reports 'ATOM record too short'"
else
    bad "A: stderr missing 'ATOM record too short' (got: ${a_err})"
fi
if printf '%s' "${a_err}" | grep -q "bad_short_atom\.mol2:${short_line}:"; then
    ok "A: stderr names fixture and line ${short_line}"
else
    bad "A: expected bad_short_atom.mol2:${short_line} in stderr (got: ${a_err})"
fi
if printf '%s' "${a_err}" | grep -Eq "out_of_range|terminate"; then
    bad "A: raised an exception instead of a diagnostic"
else
    ok "A: no out_of_range / terminate"
fi
assert_eq "A: no data rows on stdout" "$(data_rows "${a_out}")" "0"

echo "== CLI: out-of-range BOND index is rejected (B23) =="
bond_line="$(awk '
    /^@<TRIPOS>ATOM/ { inatom=1; next }
    /^@<TRIPOS>BOND/ { inatom=0; inbond=1; next }
    /^@<TRIPOS>/     { inbond=0 }
    inatom && $0 !~ /^[ \t]*$/ { natoms++ }
    inbond && NF >= 4 && ($2+0 > natoms || $3+0 > natoms) { print NR; exit }
' "${FIX}/bad_bond_index.mol2")"
run_ls_capture b_out b_err "${FIX}/bad_bond_index.mol2" "${FIX}/self.mol2" -rf 0
if printf '%s' "${b_err}" | grep -q "BOND record refers to atom index out of range"; then
    ok "B: stderr reports 'BOND record refers to atom index out of range'"
else
    bad "B: stderr missing BOND-index diagnostic (got: ${b_err})"
fi
if printf '%s' "${b_err}" | grep -q "bad_bond_index\.mol2:${bond_line}:"; then
    ok "B: stderr names fixture and line ${bond_line}"
else
    bad "B: expected bad_bond_index.mol2:${bond_line} in stderr (got: ${b_err})"
fi
if printf '%s' "${b_err}" | grep -q "valid 1-"; then
    ok "B: stderr states the valid index range"
else
    bad "B: stderr missing valid range (got: ${b_err})"
fi
assert_eq "B: no data rows on stdout" "$(data_rows "${b_out}")" "0"

# ---------------------------------------------------------------------------
# C: column values for a zero-atom ligand (B24).
# Defined zeros: PC-score8Q/8T and JaccardR.  NaN: both Pval columns (log(0))
# and RMSD_lb (the -1.0 sentinel is printed as NaN).
# ---------------------------------------------------------------------------
echo "== CLI: zero-atom ligand column values (B24) =="
c_out_f="$(mktemp)"; c_err_f="$(mktemp)"
"${BIN}" "${FIX}/empty_ligand.mol2" "${FIX}/self.mol2" -rf 0 >"${c_out_f}" 2>"${c_err_f}"
c_rc=$?
c_row="$(awk '$1=="EMPTYLIG" {print $3, $4, $5, $6, $7, $8; exit}' "${c_out_f}")"
c_err="$(cat "${c_err_f}")"
rm -f "${c_out_f}" "${c_err_f}"

if [ "${c_rc}" -eq 139 ] || [ "${c_rc}" -eq 127 ]; then
    bad "C: crashed (exit ${c_rc})"
else
    ok "C: no crash (exit ${c_rc})"
fi
assert_eq "C: PC-score8Q" "$(printf '%s
' "${c_row}" | awk '{print $1}')" "0.0000"
assert_eq "C: PC-score8T" "$(printf '%s
' "${c_row}" | awk '{print $2}')" "0.0000"
assert_eq "C: Pval_8Q" "$(printf '%s
' "${c_row}" | awk '{print $3}')" "nan"
assert_eq "C: Pval_8T" "$(printf '%s
' "${c_row}" | awk '{print $4}')" "nan"
assert_eq "C: JaccardR" "$(printf '%s
' "${c_row}" | awk '{print $5}')" "0.0000"
assert_eq "C: RMSD_lb" "$(printf '%s
' "${c_row}" | awk '{print $6}')" "nan"
if printf '%s' "${c_err}" | grep -q "has no atoms"; then
    ok "C: stderr reports 'has no atoms'"
else
    bad "C: stderr missing 'has no atoms' (got: ${c_err})"
fi

finish
