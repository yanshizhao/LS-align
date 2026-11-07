#!/usr/bin/env bash
# Black-box: non-finite mol2 input is rejected before it can reach the
# scoring chain (F2).
#   N: a NaN ATOM coordinate -> stderr "non-finite coordinate", nonzero exit
#   O: a bond whose squared length overflows to Inf -> "non-finite bond length"
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

data_rows() {   # $1 = stdout text
    printf '%s\n' "$1" | awk 'NR>2 && $3 ~ /^[0-9]+\.[0-9]+$/ {n++} END{print n+0}'
}

echo "== CLI: NaN coordinate is rejected (F2) =="
run_ls_capture n_out n_err "${FIX}/nan_coord.mol2" "${FIX}/self.mol2" -rf 0
if [ "${RC}" -ne 0 ]; then
    ok "f2_nan_coord_exits: exit code ${RC} (nonzero)"
else
    bad "f2_nan_coord_exits: exit 0 (expected nonzero)"
fi
if printf '%s' "${n_err}" | grep -q "non-finite coordinate"; then
    ok "f2_nan_coord_exits: stderr reports 'non-finite coordinate'"
else
    bad "f2_nan_coord_exits: stderr missing 'non-finite coordinate' (got: ${n_err})"
fi
if printf '%s' "${n_err}" | grep -q "nan_coord\.mol2:"; then
    ok "f2_nan_coord_exits: stderr names the fixture"
else
    bad "f2_nan_coord_exits: stderr does not name nan_coord.mol2 (got: ${n_err})"
fi
assert_eq "f2_nan_coord_exits: no data rows on stdout" "$(data_rows "${n_out}")" "0"

echo "== CLI: non-finite bond length is rejected (F2) =="
run_ls_capture o_out o_err "${FIX}/nan_bond.mol2" "${FIX}/self.mol2" -rf 0
if [ "${RC}" -ne 0 ]; then
    ok "f2_nan_bond_exits: exit code ${RC} (nonzero)"
else
    bad "f2_nan_bond_exits: exit 0 (expected nonzero)"
fi
if printf '%s' "${o_err}" | grep -q "non-finite bond length"; then
    ok "f2_nan_bond_exits: stderr reports 'non-finite bond length'"
else
    bad "f2_nan_bond_exits: stderr missing 'non-finite bond length' (got: ${o_err})"
fi
assert_eq "f2_nan_bond_exits: no data rows on stdout" "$(data_rows "${o_out}")" "0"

finish
