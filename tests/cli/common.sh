#!/usr/bin/env bash
# common.sh -- shared helpers for the LS-align black-box CLI tests.
#
# IMPORTANT: the LS-align binary ALWAYS exits with status 1 (LSalign.cpp:453),
# even on success.  So we never test $? ; we parse stdout instead.

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${LSALIGN_BIN:-${TESTS_DIR}/build/LSalign_under_test}"
EX="${TESTS_DIR}/../example"
FIX="${TESTS_DIR}/fixtures"

CLI_PASS=0
CLI_FAIL=0

# Run the binary (ignore its always-1 exit) and emit stdout only.
run_ls() {
    "${BIN}" "$@" 2>/dev/null
    return 0
}

# Extract the data row for a named query ligand: prints "PC8Q PC8T Pval8Q Pval8T Jaccard RMSD_lb Qsize Tsize".
row_for() {   # $1 = query ligand name
    run_ls "${@:2}" | awk -v name="$1" '$1==name {print $3, $4, $5, $6, $7, $8, $9, $10; exit}'
}

# Count data rows (lines after the 2 header lines whose 3rd column is numeric).
count_rows() {
    run_ls "$@" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}'
}

# assert_alive LABEL args...  -- assert the binary ran to completion.
# The binary always exits 1 on success (LSalign.cpp:453), so success is judged by
# "not a crash/load failure" (exit 139 = SIGSEGV, 127 = load/double-free abort)
# plus at least one data row on stdout.
assert_alive() {
    local label="$1"; shift
    local out rc n
    out="$("${BIN}" "$@" 2>/dev/null)"; rc=$?
    if [ "${rc}" -eq 139 ] || [ "${rc}" -eq 127 ]; then
        bad "${label}: crashed (exit ${rc})"
        return 1
    fi
    n="$(printf '%s\n' "${out}" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}')"
    if [ "${n}" -gt 0 ]; then
        ok "${label}"
    else
        bad "${label}: no data rows"
        return 1
    fi
}

# run_ls_capture OUT_VAR ERR_VAR args...  -- run the binary capturing stdout and
# stderr separately (the trailing newline of each stream is stripped).  Always
# returns 0, like run_ls.
run_ls_capture() {
    local out_var="$1"; local err_var="$2"; shift 2
    local out_f err_f
    out_f="$(mktemp)"; err_f="$(mktemp)"
    "${BIN}" "$@" >"${out_f}" 2>"${err_f}"
    eval "${out_var}=\$(cat \"${out_f}\")"
    eval "${err_var}=\$(cat \"${err_f}\")"
    rm -f "${out_f}" "${err_f}"
    return 0
}

ok()   { CLI_PASS=$((CLI_PASS+1)); echo "[  OK  ] $1"; }
bad()  { CLI_FAIL=$((CLI_FAIL+1)); echo "[ FAIL] $1"; }

# assert_close LABEL ACTUAL EXPECTED TOL
assert_close() {
    awk -v a="$2" -v e="$3" -v t="$4" -v label="$1" 'BEGIN{
        d = (a>e)? a-e : e-a;
        if (d <= t) { printf "[  OK  ] %s (=%.4f)\n", label, a; exit 0 }
        else        { printf "[ FAIL] %s (got %.4f, want %.4f +-%.4f)\n", label, a, e, t; exit 1 }
    }'
    if [ $? -eq 0 ]; then CLI_PASS=$((CLI_PASS+1)); else CLI_FAIL=$((CLI_FAIL+1)); fi
}

# assert_inrange LABEL ACTUAL LO HI
assert_inrange() {
    awk -v a="$2" -v lo="$3" -v hi="$4" -v label="$1" 'BEGIN{
        if (a>=lo && a<=hi) { printf "[  OK  ] %s (=%.4f in [%.3f,%.3f])\n", label, a, lo, hi; exit 0 }
        else                { printf "[ FAIL] %s (=%.4f NOT in [%.3f,%.3f])\n", label, a, lo, hi; exit 1 }
    }'
    if [ $? -eq 0 ]; then CLI_PASS=$((CLI_PASS+1)); else CLI_FAIL=$((CLI_FAIL+1)); fi
}

# assert_eq LABEL ACTUAL EXPECTED   (integer/string equality)
assert_eq() {
    if [ "$2" = "$3" ]; then CLI_PASS=$((CLI_PASS+1)); echo "[  OK  ] $1 (= $2)";
    else CLI_FAIL=$((CLI_FAIL+1)); echo "[ FAIL] $1 (got '$2', want '$3')"; fi
}

finish() {
    echo "----- CLI summary: ${CLI_PASS} passed, ${CLI_FAIL} failed -----"
    [ "${CLI_FAIL}" -eq 0 ]
}
