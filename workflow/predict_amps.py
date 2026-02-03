#!/usr/bin/env python3
"""
predict_amps.py
Run Macrel on a peptide FASTA and write a normalized TSV of AMP predictions.

Macrel peptide mode:
  macrel peptides --fasta <input.faa> --output <outdir>
Outputs include a "prediction" table with AMP probability + hemolysis prediction. :contentReference[oaicite:1]{index=1}
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def run(cmd: list[str], cwd: Optional[Path] = None) -> None:
    print(f"[predict_amps.py] $ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def have_exe(name: str) -> bool:
    return shutil.which(name) is not None


def macrel_supports_threads() -> bool:
    # We inspect help output to see if --threads exists (robust across versions)
    try:
        p = subprocess.run(
            ["macrel", "peptides", "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=True,
        )
        return "--threads" in p.stdout
    except Exception:
        return False


def find_prediction_file(outdir: Path) -> Path:
    """
    Macrel writes outputs into the output folder. The docs refer to an
    'expected.prediction' file in tests. We locate any file containing 'prediction'
    in its name, preferring a non-empty tabular file. :contentReference[oaicite:2]{index=2}
    """
    if not outdir.exists():
        die(f"Macrel output dir not found: {outdir}")

    candidates = []
    for p in outdir.rglob("*"):
        if p.is_file() and "prediction" in p.name.lower():
            # skip logs if any
            if p.suffix.lower() in {".log", ".err", ".out"}:
                continue
            if p.stat().st_size == 0:
                continue
            candidates.append(p)

    if not candidates:
        # fallback: sometimes tools write "predictions" etc.
        for p in outdir.rglob("*"):
            if p.is_file() and re.search(r"predic", p.name.lower()):
                if p.stat().st_size > 0:
                    candidates.append(p)

    if not candidates:
        die(f"Could not find any prediction table inside: {outdir}")

    # prefer smallest depth + most likely table extensions
    preferred_ext = [".tsv", ".txt", ".csv", ".prediction", ""]
    candidates.sort(
        key=lambda x: (
            0 if x.suffix.lower() in preferred_ext else 1,
            len(x.parts),
            x.stat().st_size,
        )
    )
    return candidates[0]


def sniff_sep(path: Path) -> str:
    # Try to detect delimiter
    sample = path.read_text(errors="replace").splitlines()[:20]
    sample_text = "\n".join(sample)
    try:
        dialect = csv.Sniffer().sniff(sample_text, delimiters="\t,;")
        return dialect.delimiter
    except Exception:
        # Macrel prediction tables are typically tab-separated
        return "\t"


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Macrel docs state the main output has 6 columns:
      accession, peptide sequence, composition/structure classification,
      AMP probability, hemolysis prediction, hemolysis probability :contentReference[oaicite:3]{index=3}
    We normalize into stable names:
      peptide_id, peptide_seq, amp_prob, amp_pred, hemo_pred, hemo_prob, (optional) class fields
    """
    # Make lowercase & strip
    cols = [str(c).strip() for c in df.columns]
    df.columns = cols

    # If file has no header (rare), create one by position
    if all(re.fullmatch(r"\d+", c) for c in cols) or len(cols) == 6 and cols[0] == 0:
        # not likely; but keep safe
        pass

    # Try to map known/likely names
    colmap = {}
    for c in df.columns:
        lc = c.lower()
        if lc in {"accession", "access_code", "access code", "id", "name", "seqid"}:
            colmap[c] = "peptide_id"
        elif lc in {"sequence", "peptide", "peptide_sequence", "peptide sequence", "seq"}:
            colmap[c] = "peptide_seq"
        elif "prob" in lc and ("amp" in lc or "antimicrobial" in lc):
            colmap[c] = "amp_prob"
        elif ("amp" in lc or "antimicrobial" in lc) and "prob" not in lc:
            # could be AMP class label
            colmap[c] = "amp_pred"
        elif ("hemo" in lc or "hemol" in lc) and "prob" in lc:
            colmap[c] = "hemo_prob"
        elif ("hemo" in lc or "hemol" in lc) and "prob" not in lc:
            colmap[c] = "hemo_pred"
        elif "cation" in lc or "anion" in lc:
            colmap[c] = "charge_class"
        elif "disulf" in lc or "linear" in lc or "structure" in lc:
            colmap[c] = "structure_class"
        elif "class" in lc or "composition" in lc:
            colmap[c] = "macrel_class"

    df = df.rename(columns=colmap)

    # If still missing, fall back to positional assumptions (6-col schema)
    needed = {"peptide_id", "peptide_seq", "amp_prob", "hemo_pred", "hemo_prob"}
    if not needed.issubset(df.columns):
        if df.shape[1] >= 6:
            # best-effort positional mapping
            # 0: accession, 1: seq, 2: class, 3: amp_prob, 4: hemo_pred, 5: hemo_prob
            tmp = df.copy()
            tmp = tmp.rename(
                columns={
                    df.columns[0]: "peptide_id",
                    df.columns[1]: "peptide_seq",
                    df.columns[2]: "macrel_class",
                    df.columns[3]: "amp_prob",
                    df.columns[4]: "hemo_pred",
                    df.columns[5]: "hemo_prob",
                }
            )
            df = tmp

    # Ensure core fields exist
    if "peptide_id" not in df.columns:
        die("Could not identify peptide_id column in Macrel output.")
    if "peptide_seq" not in df.columns:
        die("Could not identify peptide_seq column in Macrel output.")
    if "amp_prob" not in df.columns:
        die("Could not identify amp_prob column in Macrel output.")
    if "hemo_pred" not in df.columns:
        # not fatal; some versions may omit hemolysis
        df["hemo_pred"] = pd.NA
    if "hemo_prob" not in df.columns:
        df["hemo_prob"] = pd.NA

    # AMP prediction column:
    # Macrel output table (default) contains only predicted AMPs (p>0.5). :contentReference[oaicite:4]{index=4}
    if "amp_pred" not in df.columns:
        df["amp_pred"] = True

    # Clean types
    df["amp_prob"] = pd.to_numeric(df["amp_prob"], errors="coerce")
    df["hemo_prob"] = pd.to_numeric(df["hemo_prob"], errors="coerce")

    # Normalize hemo_pred to bool-ish if possible
    if df["hemo_pred"].dtype == object:
        df["hemo_pred"] = df["hemo_pred"].astype(str).str.lower().map(
            {"true": True, "false": False, "1": True, "0": False, "yes": True, "no": False}
        )

    # Final columns order
    keep = [
        "peptide_id",
        "peptide_seq",
        "amp_pred",
        "amp_prob",
        "hemo_pred",
        "hemo_prob",
    ]
    # add any optional classification fields
    for extra in ["macrel_class", "charge_class", "structure_class"]:
        if extra in df.columns:
            keep.append(extra)

    return df[keep]


def main() -> None:
    ap = argparse.ArgumentParser(description="Predict AMPs using Macrel (peptides mode).")
    ap.add_argument("fasta", type=Path, help="Input peptide FASTA (amino acids).")
    ap.add_argument("--out", type=Path, required=True, help="Output TSV path (normalized predictions).")
    ap.add_argument("--workdir", type=Path, default=None, help="Macrel output directory (default: alongside --out).")
    ap.add_argument("--threads", type=int, default=4, help="Threads for Macrel (if supported).")
    ap.add_argument("--keep-negatives", action="store_true", help="Ask Macrel to output all sequences (if supported).")
    ap.add_argument("--macrel-bin", default="macrel", help="Macrel executable name/path (default: macrel).")
    args = ap.parse_args()

    if not args.fasta.exists():
        die(f"FASTA not found: {args.fasta}")
    if not have_exe(args.macrel_bin):
        die("macrel not found on PATH. Ensure downstream env installs bioconda::macrel.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    workdir = args.workdir if args.workdir else (args.out.parent / "macrel_out_peptides")

    # If output dir exists, Macrel exits. Remove it to allow reruns.
    if workdir.exists():
        # handle symlink safely
        if workdir.is_symlink():
            workdir.unlink()
        else:
            shutil.rmtree(workdir)

    # sanity check: ensure it's really gone
    if workdir.exists():
        die(f"Could not remove existing Macrel output dir: {workdir}")

    # Build macrel command
    cmd = [args.macrel_bin, "peptides", "--fasta", str(args.fasta), "--output", str(workdir)]
    # Optional flags if supported
    if args.keep_negatives:
        cmd.append("--keep-negatives")
    if macrel_supports_threads():
        cmd += ["--threads", str(args.threads)]

    # IMPORTANT: Macrel fails if output directory exists.
    # Ensure it does NOT exist right before running.
    if workdir.exists():
        if workdir.is_symlink():
            workdir.unlink()
        else:
            shutil.rmtree(workdir)

    run(cmd)


    pred_file = find_prediction_file(workdir)
    sep = sniff_sep(pred_file)
    df = pd.read_csv(pred_file, sep=sep, engine="python")

    norm = normalize_columns(df)

    # Add provenance
    try:
        ver = subprocess.run([args.macrel_bin, "--version"], stdout=subprocess.PIPE, text=True, check=True).stdout.strip()
    except Exception:
        ver = "unknown"
    norm.insert(0, "macrel_version", ver)
    norm.insert(1, "source_fasta", str(args.fasta))

    norm.to_csv(args.out, sep="\t", index=False)
    print(f"[predict_amps.py] Wrote: {args.out}", file=sys.stderr)
    print(f"[predict_amps.py] Macrel table: {pred_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
