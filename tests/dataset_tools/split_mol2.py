#!/usr/bin/env python3
"""
split_mol2.py -- split a multi-ligand DUD-E mol2 file into single-ligand mol2
files, or keep only the first K ligand blocks.

Molecules are delimited by the line "@<TRIPOS>MOLECULE".  The LS-align binary
requires exactly one ligand per query/template file, so multi-ligand DUD-E
"actives_without_duplicate.mol2" / "decoys_without_duplicate.mol2" must be split.

Usage:
    python split_mol2.py INPUT.mol2 OUTDIR [--limit K] [--name-index I]

Options:
    --limit K        keep at most K ligand blocks (default: all)
    --name-index I   only emit the I-th (0-based) ligand block (single output)
"""
import argparse
import os
import re
import sys

HEADER = "@<TRIPOS>MOLECULE"


def split_blocks(text):
    """Return a list of (ligname, block_text) tuples."""
    lines = text.splitlines()
    blocks = []
    cur = []
    for ln in lines:
        if ln.startswith(HEADER):
            if cur:
                blocks.append(cur)
            cur = [ln]
        else:
            if cur:
                cur.append(ln)
    if cur:
        blocks.append(cur)

    out = []
    for blk in blocks:
        # molecule name is the line immediately after the MOLECULE header
        name = ""
        for i, ln in enumerate(blk):
            if ln.startswith(HEADER) and i + 1 < len(blk):
                name = blk[i + 1].strip()
                break
        out.append((name, "\n".join(blk) + "\n"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("outdir")
    ap.add_argument("--limit", type=int, default=0, help="keep at most K blocks (0 = all)")
    ap.add_argument("--name-index", type=int, default=-1, help="emit only the I-th block (0-based)")
    args = ap.parse_args()

    with open(args.input, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    blocks = split_blocks(text)
    os.makedirs(args.outdir, exist_ok=True)

    if args.name_index >= 0:
        if args.name_index >= len(blocks):
            sys.stderr.write("index %d out of range (%d blocks)\n" % (args.name_index, len(blocks)))
            return 1
        blocks = [blocks[args.name_index]]

    if args.limit and len(blocks) > args.limit:
        blocks = blocks[:args.limit]

    written = 0
    for i, (name, blk) in enumerate(blocks):
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name) or ("lig%d" % i)
        path = os.path.join(args.outdir, "%03d_%s.mol2" % (i, safe))
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(blk)
        written += 1

    sys.stderr.write("split %d ligand blocks -> %s\n" % (written, args.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
