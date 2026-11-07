# LS-align test suite

Read-only test harness for the LS-align source (`../src`). **No source files are
modified.** See `D:/qlab/LSalign_src/源码缺陷与测试报告.md` for the defect catalog
and full design.

## Quick start

```bash
make test          # unit + cli  (default CI target; no DUD-E dataset, no ASan)
make unit          # white-box unit tests (104 passed, 0 xfail, 0 xpass)
make cli           # black-box CLI tests (8 scripts)
make golden        # regenerate golden snapshots
make check-golden  # diff current output vs golden
make output-diff   # output_diff: stdout/stderr compared separately, 11 cases
```

Optional (needs the ~6.9 GB DUD-E dataset at `../dataset`, referenced not copied):
```bash
bash dataset_tools/audit_ace.sh ace 200
```

## Environment notes (this machine)

- Dynamic `libstdc++` linking is **broken** here: always pass
  `-static-libstdc++ -static-libgcc` (the Makefile does).
- `make` runs recipes with an empty `TMP`/`TEMP`/`TMPDIR`; the Makefile exports
  a writable temp dir so g++ does not fall back to `C:\Windows\`.
- ASan/UBSan are unavailable (no MinGW libasan runtime). Leak checks use the
  portable `leak_probe.h` counter instead.
- Exit code: **0 on success, 1 on failure** (B28).  Numeric assertions still parse
  stdout because the table is the contract; `cli/test_cli_args.sh` covers exit
  codes, invalid option values and stdout purity explicitly.

## xfail convention

`XFAIL(name)` encodes *correct* behaviour for a defect that is still open: the
CHECK fails today -> XFAIL (green), and once the defect is fixed the case turns
XPASS (red), forcing the maintainer to promote it to `TEST`.

The batch-1 defects are all fixed, so **no `XFAIL` cases remain**: the formerly
expected-failure cases in `unit/test_xfail.cpp` and `unit/test_delete_mismatch.cpp`
are now plain `TEST` regression locks. An XPASS is still a hard failure (it means
a lock was left un-promoted); a `TEST` failing means the defect has regressed.
New defect cases should be written directly as `TEST`.

## Unreachable paths (R5)

R5 removed the dead no-arg `LigRotamerDetector::generateRotamers()` (N=5, no
aligned-bond filter) and the three unused `FlexibleLSalign` constructors, so
only the six-argument constructor and the with-arg `generateRotamers(RigidLSalign*)`
(N=3) remain.  Consequence: the `> 5000` deterministic-decimation branch in
`generateRotamers(RigidLSalign*)` is **currently unreachable** from any public
caller (11^3 = 1331 max).  `unit/test_rotamer.cpp` keeps `x1_deterministic_decimation`,
which locks the public path's reproducibility and exercises the decimation
arithmetic through a test-only extended angle set (20 angles -> 8000 records);
it does not claim the branch is reachable in production.  `NRS_LEGACY_UNPARAM_BOND_NUM`
(=5) is retained only because the sampler tests still pass it explicitly.
