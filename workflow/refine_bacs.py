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


def parse_gff_cds_index(gff_path: str):
    """
    Build:
      - cds_by_contig: contig -> list of CDS dicts (for debug / later steps)
      - cds_index: normalized_feature_id -> (normalized_contig_id, start, end)
    """
    cds_by_contig = defaultdict(list)
    cds_index = {}

    with open(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
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
                    cds_id = attr.replace("ID=", "").strip()
                    break

            # Keep the full dict per contig (for later steps)
            cds_by_contig[seqid].append({"start": start_i, "end": end_i, "strand": strand, "id": cds_id})

            # Index by normalized CDS id (this is what we'll match feature_id against)
            if cds_id:
                cds_index[normalize_feature_id(cds_id)] = (normalize_contig_id(seqid), start_i, end_i)

    return cds_by_contig, cds_index


def ensure_cols(df: pd.DataFrame, cols: list[str]) -> None:
    for c in cols:
        if c not in df.columns:
            df[c] = ""


def normalize_contig_id(x: str) -> str:
    """
    Your contigs are like 'contig_1', 'contig_10', etc.
    We normalize by stripping a leading 'contig_' if present.
    """
    x = (x or "").strip()
    if x.startswith("contig_"):
        return x[len("contig_"):]
    return x


def normalize_feature_id(fid: str) -> str:
    """
    Your TSV feature_id can be like:
      'contig_1_1'
    Your GFF CDS ID is often like:
      '1_1'
    We normalize by stripping a leading 'contig_' if present.
    """
    fid = (fid or "").strip()
    if fid.startswith("contig_"):
        return fid[len("contig_"):]
    return fid


def tiara_is_bacteria(label) -> bool:
    if label is None:
        return False
    # pandas NaN shows up as float; also protect weird types
    if isinstance(label, float):
        return False
    s = str(label).strip().lower()
    if s in ("", "nan", "none"):
        return False
    # accept all prok buckets (depending on how tiara wrote it)
    return any(x in s for x in ("bacteria", "archaea", "prokarya", "prokaryote", "prokaryota"))


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
    ap.add_argument("--edge-bp", type=int, default=50,
                    help="Flag smORFs within this many bp of contig ends (default: 50).")

    # Behavior toggles
    ap.add_argument("--filter-tiara-bacteria", action="store_true", default=True,
                    help="Filter TSV rows to tiara_label in {bacteria, prokarya, archaea} (default: ON).")
    ap.add_argument("--no-filter-tiara-bacteria", dest="filter_tiara_bacteria", action="store_false",
                    help="Disable tiara_label filtering.")

    args = ap.parse_args()

    print("=======================================")
    print(f"Sample: {args.sample}")
    print("Building initial dictionaries...")
    print("=======================================")

    # Dict 1: contig lengths (keys are fasta headers, likely 'contig_1', ...)
    contig_len_raw = build_contig_len_dict(args.bac_contigs)
    print(f"[INFO] Contigs loaded: {len(contig_len_raw)}")

    for i, (k, v) in enumerate(contig_len_raw.items()):
        print(f"  {k} -> length={v}")
        if i == 4:
            break

    # ALSO make a normalized contig len dict keyed by stripped 'contig_' form ('1', '10', ...)
    contig_len_norm = {normalize_contig_id(k): v for k, v in contig_len_raw.items()}

    # Dict 2: CDS dict + CDS index
    cds_by_contig, cds_index = parse_gff_cds_index(args.prodigal_gff)
    print(f"[INFO] Contigs with CDS: {len(cds_by_contig)}")
    print(f"[INFO] CDS indexed by ID: {len(cds_index)}")

    for i, (contig, cds_list) in enumerate(cds_by_contig.items()):
        print(f"\n  {contig}: {len(cds_list)} CDS entries")
        for cds in cds_list[:3]:
            print(f"    {cds}")
        if i == 2:
            break

    print("\n[OK] Initial dicts built successfully.")
    print("Loading TSV and applying bacteria filter + coordinate mapping...")

    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str).copy()

    # Required columns in TSV
    required = ["feature_id", "contig_id", "tiara_label"]
    for r in required:
        if r not in df.columns:
            raise SystemExit(f"ERROR: required column '{r}' not found in TSV. Columns={list(df.columns)}")

    # Filter to bacteria (tiara_label)
    n0 = len(df)
    if args.filter_tiara_bacteria:
        lab = df["tiara_label"].fillna("").astype(str).str.lower().str.strip()
        df = df[lab.str.contains(r"(bacteria|archaea|prokarya|prokaryot)", regex=True)].copy()
        print(f"[INFO] Filter tiara_label=bacteria/prokarya/archaea: {n0} -> {len(df)} rows")

    # Ensure output schema exists
    ensure_cols(
        df,
        [
            "start", "end",  # we will populate these (even if they didn't exist before)
            "flag_edge", "dist_left", "dist_right",
            "flag_embedded", "flag_overlap_fraction", "host_cds_id",
            "cluster_id", "cluster_size", "cluster_env_count",
            "flag_cluster_recurrent", "flag_cross_environment",
            "flag_too_short",
            "confidence_tier"
        ],
    )

    # Map feature_id -> (contig, start, end) from Prodigal GFF
    mapped = 0
    unmapped = 0

    start_vals = []
    end_vals = []
    contig_vals = []  # normalize contig_id for downstream computations

    for _, row in df.iterrows():
        fid = row.get("feature_id", "")
        contig_id = row.get("contig_id", "")

        fid_norm = normalize_feature_id(fid)
        contig_norm = normalize_contig_id(contig_id)

        hit = cds_index.get(fid_norm)
        if hit is None:
            unmapped += 1
            # keep empty start/end, but still keep contig norm for step 1 (maybe)
            start_vals.append("")
            end_vals.append("")
            contig_vals.append(contig_norm)
        else:
            hit_contig_norm, s_i, e_i = hit
            mapped += 1
            # prefer GFF contig if available (more trustworthy)
            contig_vals.append(hit_contig_norm or contig_norm)
            start_vals.append(str(s_i))
            end_vals.append(str(e_i))

    df["contig_id"] = contig_vals
    df["start"] = start_vals
    df["end"] = end_vals

    print(f"[INFO] feature_id→GFF coord mapping: mapped={mapped}, unmapped={unmapped}")

    # Step 1: contig-edge artifact flagging (now we have start/end for mapped rows)
    print("Applying Step 1: contig-edge artifact flagging...")
    missing_contigs = 0
    bad_coords = 0
    flagged = 0

    dist_left_vals = []
    dist_right_vals = []
    flag_edge_vals = []

    for _, row in df.iterrows():
        contig_norm = (row.get("contig_id", "") or "").strip()
        s = (row.get("start", "") or "").strip()
        e = (row.get("end", "") or "").strip()

        try:
            s_i = int(float(s))
            e_i = int(float(e))
        except Exception:
            bad_coords += 1
            dist_left_vals.append("")
            dist_right_vals.append("")
            flag_edge_vals.append("")
            continue

        L = contig_len_norm.get(contig_norm)
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

    df["confidence_tier"] = df["confidence_tier"].replace("", "UNASSESSED")

    print(f"[INFO] Rows: {len(df)}")
    print(f"[INFO] Missing contigs in bac_contigs.fasta (after normalization): {missing_contigs}")
    print(f"[INFO] Rows with non-integer coords: {bad_coords}")
    print(f"[INFO] flag_edge=1: {flagged}")

    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
