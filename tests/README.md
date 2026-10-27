# LS-align test suite

Read-only test harness for the LS-align source (`../src`). **No source files are
modified.** See `D:/qlab/LSalign_src/源码缺陷与测试报告.md` for the defect catalog
and full design.

## Quick start

```bash
make test          # unit + cli  (default CI target; no DUD-E dataset, no ASan)
make unit          # white-box unit tests (46 passed, 0 xfail, 0 xpass)
make cli           # black-box CLI tests
make golden        # regenerate golden snapshots
make check-golden  # diff current output vs golden
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
- The LS-align binary **always exits 1**, so tests parse stdout, never `$?`.

## xfail convention

`XFAIL(name)` encodes *correct* behaviour for a defect that is still open: the
CHECK fails today -> XFAIL (green), and once the defect is fixed the case turns
XPASS (red), forcing the maintainer to promote it to `TEST`.

The batch-1 defects are all fixed, so **no `XFAIL` cases remain**: the formerly
expected-failure cases in `unit/test_xfail.cpp` and `unit/test_delete_mismatch.cpp`
are now plain `TEST` regression locks. An XPASS is still a hard failure (it means
a lock was left un-promoted); a `TEST` failing means the defect has regressed.
New defect cases should be written directly as `TEST`.
