#!/usr/bin/env python3
"""output_diff -- output-level byte/diff regression for LS-align.

Design (adapted from usalign-refactor-tests-framework/cLanguage2Cplus):
  baseline/ : expected outputs, committed, immutable (generate once).
  current/  : outputs of the current build (regenerated each run).
  diffs/    : unified diffs for any mismatch.

Two case modes (see cases.txt):
  exact          -> business data MUST be byte-identical (A-batch / F-series / refactors).
  expected:KEY   -> output may change, but it MUST differ from baseline and the diff
                    is left in diffs/ for human review against expected/KEY.json
                    (e.g. a constant->coverage mapping in R8).

Normalization before compare:
  * strip volatile timing lines  "Using ...-s", "Totally Using ...-s"
  * normalize CRLF / CR / repeated newlines to a single '\n'; outputs are
    written with newline="\n" so the committed baseline is LF on every OS
  (the outputs are textual; no binary cases are registered)

Usage:
  python output_diff.py --baseline      # (re)generate baseline/ from current build
  python output_diff.py                 # run & compare current vs baseline
  python output_diff.py --case NAME     # only one case
"""
import difflib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TESTS = HERE.parent
EXAMPLES = TESTS.parent / "example"
BIN = TESTS / "build" / ("LSalign_under_test.exe" if sys.platform.startswith("win") else "LSalign_under_test")
BASELINE = HERE / "baseline"
CURRENT = HERE / "current"
DIFFS = HERE / "diffs"
EXPECTED = HERE / "expected"

# Volatile timing lines: "Using 0.003-s to load mol2 files.",
# "Using 0.001-s to do LS-align.", "Totally Using 0.004-s." -- the number
# changes every run.
TIMING_RE = re.compile(r'^.*\b(Using|Totally Using)\b.*-s.*$', re.MULTILINE)


def cleanText(text):
    """Remove run-to-run volatile lines and normalize line endings."""
    text = TIMING_RE.sub('', text)                      # strip volatile timing lines
    text = re.sub(r'[\r\n]+', '\n', text)               # normalize newlines
    return text


def writeText(path, text):
    """Write UTF-8 with LF endings on every platform, so committed baselines
    stay byte-identical between Windows and Linux runs (tests/** is LF)."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def resolveWorkdir(rel):
    """Resolve a case workdir: 'fixtures' etc. live under tests/; 'example' is
    a sibling of tests/."""
    candTests = (TESTS / rel).resolve()
    if candTests.is_dir():
        return candTests
    candExample = (EXAMPLES if rel == "example" else TESTS.parent / rel).resolve()
    return candExample


def runCase(workdirRel, args, outPath):
    """Run one case and write its normalized output to outPath.

    LS-align always exits 1 (even on success), so the return code is reported
    only for diagnostics; pass/fail is decided by output comparison.
    """
    workdir = resolveWorkdir(workdirRel)
    cmd = [str(BIN)] + args.split()
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", cwd=str(workdir))
    content = cleanText(proc.stdout + proc.stderr)
    writeText(outPath, content)
    return proc.returncode


def diffText(basePath, curPath, tag):
    """Compare two normalized outputs; write diffs/<tag>.diff on mismatch."""
    oldLines = cleanText(basePath.read_text(encoding="utf-8", errors="replace")).splitlines(keepends=True)
    newLines = cleanText(curPath.read_text(encoding="utf-8", errors="replace")).splitlines(keepends=True)
    if oldLines == newLines:
        return "PASS"
    diffLines = list(difflib.unified_diff(oldLines, newLines,
                                          fromfile=str(basePath), tofile=str(curPath)))
    (DIFFS / f"{tag}.diff").write_text(''.join(diffLines), encoding="utf-8", newline="\n")
    return "FAIL"


def checkExpected(name, basePath, curPath, key):
    """For mode=expected:KEY.

    Light-weight contract: the case must (a) differ from the baseline (the
    change IS expected) and (b) expected/KEY.json must exist.  The unified diff
    is kept in diffs/<name>.diff for human review against the envelope declared
    in expected/KEY.json; this runner does not auto-validate numeric fields.
    """
    specPath = EXPECTED / f"{key}.json"
    if not specPath.exists():
        print(f"  WARNING: expected spec {specPath} missing -> treat as FAIL")
        return "FAIL"
    json.loads(specPath.read_text(encoding="utf-8"))     # reject malformed JSON
    if diffText(basePath, curPath, name) == "PASS":
        print("  FAIL: mode=expected but output identical to baseline (no change seen)")
        return "FAIL"
    return "EXPECTED-DIFF"


def loadCases():
    """Parse cases.txt: name | workdir | args (may contain spaces) | mode."""
    rows = []
    for line in (HERE / "cases.txt").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        head, mode = line.rsplit(None, 1)               # split off mode from the right
        parts = head.split(None, 2)                     # name, workdir, args
        if len(parts) < 3:
            print(f"  SKIP malformed: {line}")
            continue
        name, workdir, args = parts
        rows.append((name, workdir, args, mode))
    return rows


def main():
    onlyCase = None
    genBaseline = False
    argv = sys.argv[1:]
    i = 0
    while i < len(argv):
        if argv[i] == "--baseline":
            genBaseline = True
        elif argv[i] == "--case":
            if i + 1 >= len(argv):
                print("ERROR: --case requires a case name")
                sys.exit(2)
            onlyCase = argv[i + 1]
            i += 1
        else:
            print(f"ERROR: unknown argument: {argv[i]}")
            sys.exit(2)
        i += 1

    if not BIN.exists():
        print(f"ERROR: binary not found: {BIN}\nRun `make` in tests/ first.")
        sys.exit(1)
    for d in (CURRENT, DIFFS):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    BASELINE.mkdir(parents=True, exist_ok=True)

    rows = loadCases()
    total = passed = failed = warned = 0
    for name, workdir, args, mode in rows:
        if onlyCase and name != onlyCase:
            continue
        total += 1
        cur = CURRENT / f"{name}.out"
        rc = runCase(workdir, args, cur)
        base = BASELINE / f"{name}.out"

        if genBaseline:
            writeText(base, cur.read_text(encoding="utf-8"))
            print(f"  [BASELINE] {name} (rc={rc})")
            continue

        if not base.exists():
            print(f"  [MISSING] {name}: no baseline (run --baseline first)")
            failed += 1
            continue

        if mode == "exact":
            res = diffText(base, cur, name)
            if res == "PASS":
                print(f"  [  OK  ] {name} (rc={rc})")
                passed += 1
            else:
                print(f"  [ FAIL] {name} -> diffs/{name}.diff")
                failed += 1
        elif mode.startswith("expected:"):
            key = mode.split(":", 1)[1]
            res = checkExpected(name, base, cur, key)
            if res == "EXPECTED-DIFF":
                print(f"  [ WARN ] {name}: change matches expected:{key} (review diffs/{name}.diff)")
                warned += 1
            else:
                print(f"  [ FAIL] {name} -> diffs/{name}.diff")
                failed += 1
        else:
            print(f"  [ FAIL] {name}: unknown mode {mode}")
            failed += 1

    if genBaseline:
        print(f"\n=== baseline written: {total} cases ===")
        return
    if onlyCase and total == 0:
        print(f"ERROR: no case named '{onlyCase}' in cases.txt")
        sys.exit(2)
    print(f"\n=== output_diff: total={total} PASS={passed} WARN={warned} FAIL={failed} ===")
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
