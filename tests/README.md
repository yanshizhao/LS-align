# LS-align test suite

Read-only test harness for the LS-align source (`../src`). **No source files are
modified.** See `D:/qlab/LSalign_src/源码缺陷与测试报告.md` for the defect catalog
and full design.

## Quick start

```bash
make test          # unit + cli  (default CI target; no DUD-E dataset, no ASan)
make unit          # white-box unit tests (31 passed, 5 xfail expected)
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

`XFAIL(name)` encodes *correct* behaviour. Today the code is buggy -> the CHECK
fails -> XFAIL (green). When fixed it becomes XPASS (red), forcing you to
promote it to `TEST`. This makes regressions impossible to hide.
