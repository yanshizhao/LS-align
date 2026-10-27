#!/usr/bin/env bash
# check_golden.sh -- regenerate the current CLI output and diff against the
# committed golden snapshots.  Deterministic: rigid mode (-rf 0) only.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_DIR="$(cd "${HERE}/.." && pwd)"
BIN="${TESTS_DIR}/build/LSalign_under_test"
EX="${TESTS_DIR}/../example"
G="${TESTS_DIR}/golden"
TMPD="$(mktemp -d)"
rc=0

gen() {   # $1 = output file, rest = args
    # strip the three variable timing lines ("Using ...-s", "Totally Using ...")
    "${BIN}" "${@:2}" 2>/dev/null \
      | grep -vE '^(Using|Totally Using) ' > "$1" || true   # binary always exits 1
}

check_one() {   # $1 = golden file, $2 = label, rest = args
    local gold="$1"; local label="$2"; shift 2
    gen "${TMPD}/cur.txt" "$@"
    if diff -q "${gold}" "${TMPD}/cur.txt" >/dev/null 2>&1; then
        echo "[  OK  ] golden match: ${label}"
    else
        echo "[ FAIL] golden DIFF: ${label}"
        diff "${gold}" "${TMPD}/cur.txt" | head -20
        rc=1
    fi
}

if [ ! -f "${G}/ace_rigid.txt" ]; then
    echo "golden/ace_rigid.txt missing -- run 'make golden' first"; exit 1
fi

check_one "${G}/ace_rigid.txt"  "ace_ligands x xtal-lig (-rf 0)" "${EX}/ace_ligands.mol2" "${EX}/xtal-lig.mol2" -rf 0
if [ -f "${G}/self_copy.txt" ]; then
    check_one "${G}/self_copy.txt" "self x copy (-rf 0)" "${TESTS_DIR}/fixtures/self.mol2" "${TESTS_DIR}/fixtures/copy.mol2" -rf 0
fi

rm -rf "${TMPD}"
[ "${rc}" -eq 0 ] && echo "----- golden: all snapshots match -----" || echo "----- golden: DIFFERENCES FOUND -----"
exit ${rc}
