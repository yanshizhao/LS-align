#!/usr/bin/env bash
# Black-box: mol2 parsing robustness (F5) and PDB atom-name output (B29).
#   P: -o writes original atom names to 13-16 and the element to 77-78
#   Q: an unparseable ATOM record is rejected, not silently zeroed
#   R: an unknown extension falls back to the mol2 parser with a warning
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

TMPD="$(mktemp -d)"

echo "== CLI: -o writes atom names and element columns (B29) =="
run_ls "${EX}/xtal-lig.mol2" "${FIX}/self.mol2" -rf 0 -o "${TMPD}/sup.pdb" >/dev/null
first="$(grep -m1 '^HETATM' "${TMPD}/sup.pdb")"
assert_eq "P: first atom name (13-16)" "$(printf '%s' "${first}" | cut -c13-16)" " O1 "
assert_eq "P: first element (77-78)"   "$(printf '%s' "${first}" | cut -c77-78)" " O"
assert_eq "P: residue stays at 18-20"  "$(printf '%s' "${first}" | cut -c18-20)" "QUE"
assert_eq "P: every HETATM has a name" \
    "$(grep '^HETATM' "${TMPD}/sup.pdb" | awk 'substr($0,13,4) ~ /^[ ]*$/ {n++} END{print n+0}')" "0"

echo "== CLI: unparseable ATOM record is rejected (F5) =="
run_ls_capture q_out q_err "${FIX}/bad_atom_coords.mol2" "${FIX}/self.mol2" -rf 0
if [ "${RC}" -ne 0 ]; then
    ok "Q: exit code ${RC} (nonzero)"
else
    bad "Q: exit 0 (expected nonzero)"
fi
if printf '%s' "${q_err}" | grep -q "cannot parse ATOM record"; then
    ok "Q: stderr reports 'cannot parse ATOM record'"
else
    bad "Q: stderr missing parse diagnostic (got: ${q_err})"
fi
assert_eq "Q: no data rows on stdout" \
    "$(printf '%s\n' "${q_out}" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}')" "0"

echo "== CLI: unknown extension falls back to mol2 (F5) =="
run_ls_capture r_out r_err "${FIX}/self.pdb" "${FIX}/self.mol2" -rf 0
assert_eq "R: exit 0 (fallback parsed)" "${RC}" "0"
if printf '%s' "${r_err}" | grep -q "unsupported file extension"; then
    ok "R: stderr reports 'unsupported file extension'"
else
    bad "R: stderr missing extension warning (got: ${r_err})"
fi
assert_eq "R: query size parsed as mol2" \
    "$(printf '%s\n' "${r_out}" | awk '$1=="SELF" {print $9; exit}')" "6"

rm -rf "${TMPD}"
finish
