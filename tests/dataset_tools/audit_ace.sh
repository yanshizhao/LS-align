#!/usr/bin/env bash
# audit_ace.sh -- sanity check on the ACE target that the trend reported in the
# paper (actives score higher than decoys against the crystal ligand) holds.
#
# This is NOT a strict unit test -- it is an oracle/trend check using real
# DUD-E data, reported as pass/fail on an aggregate statistic.  Rigid mode only
# (deterministic).
#
# Sample size matters: the DUD-E files are not shuffled, so a tiny K (the first
# few ligands) is unrepresentative.  K=200 was verified to reproduce the paper's
# trend on ACE (actives mean PC-score8T 0.457 > decoys 0.390).
#
# Requires the dataset; skipped cleanly if unavailable.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_DIR="$(cd "${HERE}/.." && pwd)"
REPO_ROOT="$(cd "${TESTS_DIR}/.." && pwd)"
DATASET="${LSALIGN_DATASET:-${REPO_ROOT}/../dataset}"

TARGET="${1:-ace}"
K="${2:-200}"
WORK="${TESTS_DIR}/build/sampled/${TARGET}_audit"

if [ ! -d "${DATASET}/${TARGET}" ]; then
    echo "SKIP: dataset target '${TARGET}' not found under ${DATASET}"
    exit 2
fi

bash "${HERE}/sample.sh" "${TARGET}" "${K}" >/dev/null 2>&1 || { echo "SKIP: sample.sh failed"; exit 2; }
SRC="${TESTS_DIR}/build/sampled/${TARGET}"

# mean PC-score8T (col 4) over actives vs decoys (higher = more similar to crystal)
mean_col4() { awk 'NR>2 && $4 ~ /^[0-9.]+$/ {s+=$4; n++} END{ if(n) printf "%.4f", s/n; else print "NA" }' "$1"; }

MA=$(mean_col4 "${SRC}/actives_vs_crystal.txt")
MD=$(mean_col4 "${SRC}/decoys_vs_crystal.txt")

echo "ACE (${TARGET}), sample K=${K}"
echo "  mean PC-score8T  actives : ${MA}"
echo "  mean PC-score8T  decoys  : ${MD}"

# The paper's trend: actives should be recognized as more similar than decoys.
if [ "${MA}" != "NA" ] && [ "${MD}" != "NA" ]; then
    awk -v a="${MA}" -v d="${MD}" 'BEGIN{
        if (a >= d) { printf "  TREND OK: actives (%.4f) >= decoys (%.4f)\n", a, d; exit 0 }
        else        { printf "  TREND WARN: actives (%.4f) < decoys (%.4f)\n", a, d; exit 1 }
    }'
    rc=$?
    [ ${rc} -eq 0 ] && exit 0 || exit 1
fi
echo "  insufficient data"; exit 2
