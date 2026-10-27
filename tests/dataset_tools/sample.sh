#!/usr/bin/env bash
# sample.sh -- build a small, deterministic working set from the DUD-E dataset
# and run LS-align over it.
#
# The 6.9 GB dataset is REFERENCED, never copied:
#   LSALIGN_DATASET  (default ../dataset relative to the repo root)
#
# For each requested target this script takes the crystal ligand as the
# reference and samples the first K actives and first K decoys (split to single
# ligands).  Rigid mode (-rf 0) is used so results are reproducible; flexible
# mode under >5000 rotamers is nondeterministic (srand(clock())).
#
# Usage:
#   bash sample.sh <target> [K]
#   bash sample.sh ace 10
#
# Outputs land in ../../build/sampled/<target>/ (git-ignored working area).
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_DIR="$(cd "${HERE}/.." && pwd)"
REPO_ROOT="$(cd "${TESTS_DIR}/.." && pwd)"

DATASET="${LSALIGN_DATASET:-${REPO_ROOT}/../dataset}"
BIN="${TESTS_DIR}/build/LSalign_under_test"
WORK="${TESTS_DIR}/build/sampled"

TARGET="${1:-ace}"
K="${2:-10}"

TDIR="${DATASET}/${TARGET}"
if [ ! -d "${TDIR}" ]; then
    echo "no such target dir: ${TDIR}"; exit 1
fi
if [ ! -x "${BIN}" ]; then
    echo "binary missing: ${BIN}  (run 'make cli' first)"; exit 1
fi

OUT="${WORK}/${TARGET}"
rm -rf "${OUT}"; mkdir -p "${OUT}/actives" "${OUT}/decoys"

echo "== splitting ${TARGET}: actives (first ${K}) =="
python "${HERE}/split_mol2.py" "${TDIR}/actives_without_duplicate.mol2" "${OUT}/actives" --limit "${K}"
echo "== splitting ${TARGET}: decoys (first ${K}) =="
python "${HERE}/split_mol2.py" "${TDIR}/decoys_without_duplicate.mol2" "${OUT}/decoys" --limit "${K}"

# cat the K single-ligand files back into one multi-ligand query file
cat "${OUT}/actives/"*.mol2 > "${OUT}/actives_K.mol2"
cat "${OUT}/decoys/"*.mol2  > "${OUT}/decoys_K.mol2"

echo "== running LS-align (rigid) over actives vs crystal ligand =="
"${BIN}" "${OUT}/actives_K.mol2" "${TDIR}/crystal_ligand.mol2" -rf 0 2>/dev/null \
    | grep -vE '^(Using|Totally Using) ' > "${OUT}/actives_vs_crystal.txt" || true

echo "== running LS-align (rigid) over decoys vs crystal ligand =="
"${BIN}" "${OUT}/decoys_K.mol2" "${TDIR}/crystal_ligand.mol2" -rf 0 2>/dev/null \
    | grep -vE '^(Using|Totally Using) ' > "${OUT}/decoys_vs_crystal.txt" || true

echo "results -> ${OUT}"
echo "  actives_vs_crystal.txt : $(wc -l < "${OUT}/actives_vs_crystal.txt") lines"
echo "  decoys_vs_crystal.txt  : $(wc -l < "${OUT}/decoys_vs_crystal.txt") lines"
