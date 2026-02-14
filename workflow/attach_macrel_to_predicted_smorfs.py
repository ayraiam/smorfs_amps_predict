#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path
import pandas as pd


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def guess_col(cols: list[str], candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Attach Macrel predictions (from predict_amps.py) onto predicted_smorfs.tsv by peptide ID."
    )
    ap.add_argument("--predicted", required=True, type=Path, help="predicted_smorfs.tsv")
    ap.add_argument("--macrel", required=True, type=Path, help="Normalized Macrel TSV from workflow/predict_amps.py")
    ap.add_argument("--out", required=True, type=Path, help="Output TSV (merged)")
    ap.add_argument("--id-col", default=None, help="ID column name in predicted_smorfs.tsv (auto-detect if omitted)")
    ap.add_argument("--seq-col", default=None, help="AA sequence column name in predicted_smorfs.tsv (auto-detect if omitted)")
    ap.add_argument(
        "--require-seq-match",
        action="store_true",
        default=False,
        help="Fail if peptide sequences do not match between tables (default: OFF)",
    )

    args = ap.parse_args()

    if not args.predicted.exists():
        die(f"Missing predicted TSV: {args.predicted}")
    if not args.macrel.exists():
        die(f"Missing macrel TSV: {args.macrel}")

    pred = pd.read_csv(args.predicted, sep="\t", dtype=str)
    mac = pd.read_csv(args.macrel, sep="\t", dtype=str)

    # Macrel normalized columns produced by your predict_amps.py
    # expected: peptide_id, peptide_seq, amp_pred, amp_prob, hemo_pred, hemo_prob (+ optional class fields)
    if "peptide_id" not in mac.columns:
        die(f"Macrel TSV missing 'peptide_id'. Columns: {list(mac.columns)}")
    if "peptide_seq" not in mac.columns:
        die(f"Macrel TSV missing 'peptide_seq'. Columns: {list(mac.columns)}")

    id_col = args.id_col or guess_col(
        list(pred.columns),
        ["peptide_id", "smorf_id", "orf_id", "id", "accession", "name", "seqid"]
    )
    if not id_col:
        die(f"Could not auto-detect ID column in predicted TSV. Columns: {list(pred.columns)}")

    seq_col = args.seq_col or guess_col(
        list(pred.columns),
        ["peptide_seq", "aa_seq", "aa", "sequence", "pep_seq", "peptide"]
    )
    if not seq_col:
        die(f"Could not auto-detect AA sequence column in predicted TSV. Columns: {list(pred.columns)}")

    # Clean whitespace
    pred[id_col] = pred[id_col].astype(str).str.strip()
    pred[seq_col] = pred[seq_col].astype(str).str.strip()
    mac["peptide_id"] = mac["peptide_id"].astype(str).str.strip()
    mac["peptide_seq"] = mac["peptide_seq"].astype(str).str.strip()

    # Check duplicates
    if pred[id_col].duplicated().any():
        dups = pred.loc[pred[id_col].duplicated(), id_col].head(10).tolist()
        die(f"predicted TSV has duplicated IDs in {id_col}. Example duplicates: {dups}")

    # Drop placeholder columns from predicted table (we'll append Macrel columns instead)
    placeholders = ["amp_pred", "amp_score", "hemolytic", "toxic", "notes"]
    pred = pred.drop(columns=[c for c in placeholders if c in pred.columns], errors="ignore")

    # Macrel might output only predicted AMPs depending on version/settings.
    # We do a left join and fill missing predictions with NA.
    merged = pred.merge(
        mac,
        how="left",
        left_on=id_col,
        right_on="peptide_id",
        suffixes=("", "_macrel")
    )

    # Sequence validation for matched rows
    # Optional sequence validation for matched rows (OFF by default)
    if args.require_seq_match:
        matched = merged["peptide_id"].notna()
        if matched.any():
            seq_pred = merged.loc[matched, seq_col].astype(str).str.strip().str.rstrip("*")
            seq_mac  = merged.loc[matched, "peptide_seq"].astype(str).str.strip().str.rstrip("*")
            mism = (seq_pred != seq_mac)
            if mism.any():
                bad = merged.loc[matched].loc[mism, [id_col, seq_col, "peptide_seq"]].head(10)
                die(
                    "Sequence mismatch for some IDs between predicted TSV and Macrel TSV.\n"
                    f"First mismatches:\n{bad.to_string(index=False)}"
                )

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote merged TSV: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

