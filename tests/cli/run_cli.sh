#!/usr/bin/env bash
# run_cli.sh -- orchestrates the black-box CLI tests.
# Exits non-zero if any individual test script reports a failure.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

rc=0
for t in test_cli_self.sh test_cli_ace.sh test_cli_d0.sh test_cli_outputs.sh test_cli_badmol2.sh; do
    echo "########## ${t} ##########"
    bash "${HERE}/${t}" || rc=1
    echo
done
exit ${rc}
