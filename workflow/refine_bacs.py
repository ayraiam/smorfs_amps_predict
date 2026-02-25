#!/usr/bin/env python3
import argparse
from collections import defaultdict

import pandas as pd
from Bio import SeqIO


def build_contig_len_dict(fasta_path: str) -> dict[str, int]:
    contig_len = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        contig_len[record.id] = len(record.seq)
    return contig_len


def build_cds_dict_from_gff(gff_path: str):
    cds_by_contig = defaultdict(list)

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            if feature_type != "CDS":
                continue

            start_i = int(start)
            end_i = int(end)

            cds_id = ""
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    cds_id = attr.replace("ID=", "")
                    break

            cds_by_contig[seqid].append({"start": start_i, "end": end_i, "strand": strand, "id": cds_id})

    return cds_by_contig


def ensure_cols(df: pd.DataFrame, cols: list[str]) -> None:
    for c in cols:
        if c not in df.columns:
            df[c] = ""


def pick_col(df: pd.DataFrame, candidates: list[str], label: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise SystemExit(
        f"ERROR: could not find a '{label}' column in TSV. "
        f"Tried: {candidates}. Available columns: {list(df.columns)}"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-tsv", required=True)
    ap.add_argument("--prodigal-gff", required=True)
    ap.add_argument("--bac-contigs", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpus", type=int, default=8)

    # Step 1 config
    ap.add_argument("--edge-bp", type=int, default=50, help="Flag smORFs within this many bp of contig ends (default: 50).")

    # Optional explicit column names (recommended if your TSV schema is stable)
    ap.add_argument("--contig-col", default="", help="Contig/seqid column name in TSV (optional).")
    ap.add_argument("--start-col", default="", help="Start coordinate column name in TSV (optional).")
    ap.add_argument("--end-col", default="", help="End coordinate column name in TSV (optional).")

    args = ap.parse_args()

    print("=======================================")
    print(f"Sample: {args.sample}")
    print("Building initial dictionaries...")
    print("=======================================")

    # Dict 1: contig lengths
    contig_len = build_contig_len_dict(args.bac_contigs)
    print(f"[INFO] Contigs loaded: {len(contig_len)}")
    for i, (k, v) in enumerate(contig_len.items()):
        print(f"  {k} -> length={v}")
        if i == 4:
            break

    # Dict 2: CDS per contig (not used yet in Step 1, but ok to keep)
    cds_by_contig = build_cds_dict_from_gff(args.prodigal_gff)
    print(f"[INFO] Contigs with CDS: {len(cds_by_contig)}")
    for i, (contig, cds_list) in enumerate(cds_by_contig.items()):
        print(f"\n  {contig}: {len(cds_list)} CDS entries")
        for cds in cds_list[:3]:
            print(f"    {cds}")
        if i == 2:
            break

    print("\n[OK] Initial dicts built successfully.")
    print("Applying Step 1: contig-edge artifact flagging...")

    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str).copy()

    # Ensure output schema exists
    ensure_cols(
        df,
        [
            "flag_edge", "dist_left", "dist_right",
            "flag_embedded", "flag_overlap_fraction", "host_cds_id",
            "cluster_id", "cluster_size", "cluster_env_count",
            "flag_cluster_recurrent", "flag_cross_environment",
            "flag_too_short",
            "confidence_tier"
        ],
    )

    # Resolve columns for Step 1
    contig_candidates = ["contig", "contig_id", "seqid", "chrom", "scaffold", "contig_name"]
    start_candidates = ["start", "start_bp", "begin", "begin_bp", "orf_start", "gene_start"]
    end_candidates = ["end", "end_bp", "stop", "stop_bp", "orf_end", "gene_end"]

    contig_col = args.contig_col or pick_col(df, contig_candidates, "contig/seqid")
    start_col = args.start_col or pick_col(df, start_candidates, "start")
    end_col = args.end_col or pick_col(df, end_candidates, "end")

    print(f"[INFO] Using columns: contig={contig_col}, start={start_col}, end={end_col}, edge_bp={args.edge_bp}")

    missing_contigs = 0
    bad_coords = 0
    flagged = 0

    # Compute
    dist_left_vals = []
    dist_right_vals = []
    flag_edge_vals = []

    for _, row in df.iterrows():
        contig = row.get(contig_col, "")
        s = row.get(start_col, "")
        e = row.get(end_col, "")

        try:
            s_i = int(float(s))  # tolerate "123.0"
            e_i = int(float(e))
        except Exception:
            bad_coords += 1
            dist_left_vals.append("")
            dist_right_vals.append("")
            flag_edge_vals.append("")
            continue

        L = contig_len.get(contig)
        if L is None:
            missing_contigs += 1
            dist_left_vals.append("")
            dist_right_vals.append("")
            flag_edge_vals.append("")
            continue

        dl = s_i - 1
        dr = L - e_i
        is_edge = 1 if (dl < args.edge_bp or dr < args.edge_bp) else 0
        if is_edge:
            flagged += 1

        dist_left_vals.append(str(dl))
        dist_right_vals.append(str(dr))
        flag_edge_vals.append(str(is_edge))

    df["dist_left"] = dist_left_vals
    df["dist_right"] = dist_right_vals
    df["flag_edge"] = flag_edge_vals

    # Keep tier as UNASSESSED for now (Step 7 later)
    df["confidence_tier"] = df["confidence_tier"].replace("", "UNASSESSED")

    print(f"[INFO] Rows: {len(df)}")
    print(f"[INFO] Missing contigs in bac_contigs.fasta: {missing_contigs}")
    print(f"[INFO] Rows with non-integer coords: {bad_coords}")
    print(f"[INFO] flag_edge=1: {flagged}")

    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
