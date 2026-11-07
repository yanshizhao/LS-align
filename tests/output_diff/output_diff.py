#!/usr/bin/env python3
"""output_diff -- output-level byte/diff regression for LS-align.

Design (adapted from usalign-refactor-tests-framework/cLanguage2Cplus):
  baseline/ : expected outputs, committed, immutable (generate once).
  current/  : outputs of the current build (regenerated each run).
  diffs/    : unified diffs for any mismatch.

stdout and stderr are captured, normalized and compared SEPARATELY:
  baseline/<name>.out  vs current/<name>.out   (result table + timing)
  baseline/<name>.err  vs current/<name>.err   (diagnostics)
A case passes only when every stream matches; mismatches produce
diffs/<name>.out.diff and/or diffs/<name>.err.diff.  Keeping the streams apart
is what lets B27 (diagnostics moved to stderr) be verified: the result table
must stay byte-identical while diagnostics leave stdout.

Two case modes (see cases.txt):
  exact          -> both streams MUST be byte-identical (A-batch / F-series / refactors).
  expected:KEY   -> at least one stream may change, but the output MUST differ from
                    baseline and the diffs are left in diffs/ for human review
                    against expected/KEY.json (e.g. a constant->coverage mapping in R8).

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


def runCase(workdirRel, args, outPath, errPath):
    """Run one case and write its normalized stdout/stderr to outPath/errPath.

    Output comparison decides pass/fail, as before; the return code is reported
    only for diagnostics (it is 0 on success after B28).
    """
    workdir = resolveWorkdir(workdirRel)
    cmd = [str(BIN)] + args.split()
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", cwd=str(workdir))
    writeText(outPath, cleanText(proc.stdout))
    writeText(errPath, cleanText(proc.stderr))
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


def compareStreams(name):
    """Compare stdout and stderr independently.

    Returns (allMatch, anyDifference); diffs/<name>.out.diff and/or
    diffs/<name>.err.diff are written for whichever stream differs.
    """
    outRes = diffText(BASELINE / f"{name}.out", CURRENT / f"{name}.out", f"{name}.out")
    errRes = diffText(BASELINE / f"{name}.err", CURRENT / f"{name}.err", f"{name}.err")
    return (outRes == "PASS" and errRes == "PASS"), (outRes != "PASS" or errRes != "PASS")


def checkExpected(name, key):
    """For mode=expected:KEY.

    Light-weight contract: the case must (a) differ from the baseline in at
    least one stream (the change IS expected) and (b) expected/KEY.json must
    exist.  The unified diffs are kept in diffs/<name>.{out,err}.diff for human
    review against the envelope declared in expected/KEY.json; this runner does
    not auto-validate numeric fields.
    """
    specPath = EXPECTED / f"{key}.json"
    if not specPath.exists():
        print(f"  WARNING: expected spec {specPath} missing -> treat as FAIL")
        return "FAIL"
    json.loads(specPath.read_text(encoding="utf-8"))     # reject malformed JSON
    if not compareStreams(name)[1]:
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
        curOut = CURRENT / f"{name}.out"
        curErr = CURRENT / f"{name}.err"
        rc = runCase(workdir, args, curOut, curErr)
        baseOut = BASELINE / f"{name}.out"
        baseErr = BASELINE / f"{name}.err"

        if genBaseline:
            writeText(baseOut, curOut.read_text(encoding="utf-8"))
            writeText(baseErr, curErr.read_text(encoding="utf-8"))
            print(f"  [BASELINE] {name} (rc={rc})")
            continue

        if not baseOut.exists() or not baseErr.exists():
            print(f"  [MISSING] {name}: no baseline (run --baseline first)")
            failed += 1
            continue

        if mode == "exact":
            if compareStreams(name)[0]:
                print(f"  [  OK  ] {name} (rc={rc})")
                passed += 1
            else:
                print(f"  [ FAIL] {name} -> diffs/{name}.out.diff / {name}.err.diff")
                failed += 1
        elif mode.startswith("expected:"):
            key = mode.split(":", 1)[1]
            res = checkExpected(name, key)
            if res == "EXPECTED-DIFF":
                print(f"  [ WARN ] {name}: change matches expected:{key} (review diffs/{name}.{{out,err}}.diff)")
                warned += 1
            else:
                print(f"  [ FAIL] {name}")
                failed += 1
        else:
            print(f"  [ FAIL] {name}: unknown mode {mode}")
            failed += 1

    if genBaseline:
        print(f"\n=== baseline written: {total} cases (.out + .err) ===")
        return
    if onlyCase and total == 0:
        print(f"ERROR: no case named '{onlyCase}' in cases.txt")
        sys.exit(2)
    print(f"\n=== output_diff: total={total} PASS={passed} WARN={warned} FAIL={failed} ===")
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
