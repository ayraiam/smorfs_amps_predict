#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

def build_contig_len_dict(fasta_path: str) -> dict[str, int]:
    contig_len = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        contig_len[record.id] = len(record.seq)
    return contig_len


def parse_gff_cds_index(gff_path: str):
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

            # parse attributes into a dict
            ad = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    ad[k.strip()] = v.strip()

            cds_id = ad.get("ID", "").strip()
            cds_name = ad.get("Name", "").strip()
            cds_locus = ad.get("locus_tag", "").strip()

            seqid_norm = normalize_contig_id(seqid)
            cds_by_contig[seqid_norm].append({"start": start_i, "end": end_i, "strand": strand, "id": cds_id})

            # index multiple keys -> same coords
            for key in (cds_id, cds_name, cds_locus):
                if key:
                    cds_index[key] = (seqid_norm, start_i, end_i)

    return cds_by_contig, cds_index

def parse_smorfinder_gff_index(smorf_gff_path: str):
    smorf_by_contig = defaultdict(list)
    smorf_index = {}

    if not smorf_gff_path:
        return smorf_by_contig, smorf_index

    with open(smorf_gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields

            # DO NOT filter by feature_type; SmORFinder may not use "CDS"
            fid = ""
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    fid = attr.replace("ID=", "").strip()
                    break
            if not fid:
                continue

            s_i = int(start)
            e_i = int(end)

            seqid_norm = normalize_contig_id(seqid)
            smorf_by_contig[seqid_norm].append({"id": fid, "start": s_i, "end": e_i, "strand": strand})
            smorf_index[fid] = (seqid_norm, s_i, e_i, strand)

    return smorf_by_contig, smorf_index


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

FEATURE_RE = re.compile(r"^(contig_\d+)_(\d+)$")

def parse_feature_id(fid: str):
    """
    fid like: contig_100_20
    returns: ("contig_100", 20) or (None, None) if not parseable
    """
    fid = (fid or "").strip()
    m = FEATURE_RE.match(fid)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))

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

def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    s = max(a_start, b_start)
    e = min(a_end, b_end)
    return max(0, e - s + 1)

def compute_overlaps_all_hosts(sm_start: int, sm_end: int, cds_list: list[dict]):
    """
    Returns:
      host_ids: list[str]
      fracs: list[float]   # overlap fraction per host CDS (over smORF length)
      embedded_hosts: list[str]  # hosts that fully contain smORF
    """
    sm_len = sm_end - sm_start + 1
    if sm_len <= 0:
        return [], [], []

    host_ids = []
    fracs = []
    embedded_hosts = []

    for cds in cds_list:
        cs, ce = int(cds["start"]), int(cds["end"])
        ov = overlap_len(sm_start, sm_end, cs, ce)
        if ov <= 0:
            continue
        frac = ov / sm_len
        host_ids.append(str(cds.get("id", "")))
        fracs.append(frac)

        if cs <= sm_start and ce >= sm_end:
            embedded_hosts.append(str(cds.get("id", "")))

    return host_ids, fracs, embedded_hosts

def compute_step2_embedded_and_overlap(df: pd.DataFrame, cds_by_contig, min_host_extra_bp: int = 0):
    """
    For each smORF row in df (needs contig_id/start/end), compute:
      - flag_embedded (1/0)
      - host_cds_id (best host)
      - flag_overlap_fraction (max overlap with any CDS, as fraction of smORF length)
    """
    flag_embedded_vals = []
    host_cds_id_vals = []
    overlap_frac_vals = []

    embedded_n = 0
    checked_n = 0
    missing_coords = 0
    missing_contig = 0

    for _, row in df.iterrows():
        contig = (row.get("contig_id", "") or "").strip()
        s = (row.get("start", "") or "").strip()
        e = (row.get("end", "") or "").strip()

        try:
            s_i = int(float(s))
            e_i = int(float(e))
        except Exception:
            missing_coords += 1
            flag_embedded_vals.append("")
            host_cds_id_vals.append("")
            overlap_frac_vals.append("")
            continue

        if not contig:
            missing_contig += 1
            flag_embedded_vals.append("")
            host_cds_id_vals.append("")
            overlap_frac_vals.append("")
            continue

        cds_list = cds_by_contig.get(contig)
        if not cds_list:
            # no CDS on contig
            flag_embedded_vals.append("0")
            host_cds_id_vals.append("")
            overlap_frac_vals.append("0")
            continue

        checked_n += 1
        sm_len = e_i - s_i + 1

        # 1) embedded check: any CDS fully contains smORF AND is longer
        best_host = None
        best_host_len = None

        # 2) overlap fraction: max overlap with any CDS
        max_ov = 0

        for cds in cds_list:
            cs = int(cds["start"])
            ce = int(cds["end"])
            cds_len = ce - cs + 1

            ov = overlap_len(s_i, e_i, cs, ce)
            if ov > max_ov:
                max_ov = ov

            # containment
            if cs <= s_i and ce >= e_i:
                # require host longer by at least min_host_extra_bp
                if (cds_len - sm_len) >= min_host_extra_bp:
                    if best_host is None or cds_len < best_host_len:
                        best_host = cds
                        best_host_len = cds_len

        if best_host is not None:
            embedded_n += 1
            flag_embedded_vals.append("1")
            host_cds_id_vals.append(str(best_host.get("id", "")))
        else:
            flag_embedded_vals.append("0")
            host_cds_id_vals.append("")

        frac = (max_ov / sm_len) if sm_len > 0 else 0.0
        # keep a compact representation
        overlap_frac_vals.append(f"{frac:.3f}")

    return {
        "flag_embedded": flag_embedded_vals,
        "host_cds_id_step2": host_cds_id_vals,  # we’ll merge carefully to not clobber your existing host_cds_id unless you want it
        "flag_overlap_fraction": overlap_frac_vals,
        "stats": {
            "checked_rows": checked_n,
            "embedded_rows": embedded_n,
            "missing_coords": missing_coords,
            "missing_contig": missing_contig,
        }
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-tsv", required=True)
    ap.add_argument("--prodigal-gff", default="",
                help="Prodigal bacterial genes GFF. If missing/empty, Step2 is skipped.")
    ap.add_argument("--bac-contigs", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpus", type=int, default=8)

    # Step 1 config
    ap.add_argument("--edge-bp", type=int, default=50,
                    help="Flag smORFs within this many bp of contig ends (default: 50).")

    # Step 2 config
    ap.add_argument("--min-host-extra-bp", type=int, default=0,
                    help="Require host CDS to be at least this many bp longer than smORF to call embedded (default: 0).")

    # Behavior toggles
    ap.add_argument(
        "--no-filter-tiara-bacteria",
        dest="filter_tiara_bacteria",
        action="store_false",
        help="Disable tiara_label filtering (default: filter is ON).",
    )

    ap.add_argument(
    "--smorfinder-gff", default="",
    help="SmORFinder GFF (smorf_output.gff). Used to map smorfinder feature_id -> contig/start/end."
    )

    ap.set_defaults(filter_tiara_bacteria=True)

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

    smorf_by_contig, smorf_index = parse_smorfinder_gff_index(args.smorfinder_gff)
    print(f"[INFO] SmORFinder contigs with smORFs: {len(smorf_by_contig)}")
    print(f"[INFO] SmORFinder features indexed: {len(smorf_index)}")

    for i, (contig, sm_list) in enumerate(smorf_by_contig.items()):
        print(f"\n  {contig}: {len(sm_list)} smORFs")
        for sm in sm_list[:3]:
             print(f"    {sm}")
        if i == 2:
            break

    print("\n[OK] Initial dicts built successfully.")
    print("Loading TSV and applying bacteria filter + coordinate mapping...")

    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str).copy()

    # Required columns in TSV
    required = ["source", "feature_id", "contig_id", "tiara_label"]
    for r in required:
        if r not in df.columns:
            raise SystemExit(f"ERROR: required column '{r}' not found in TSV. Columns={list(df.columns)}")

    # Filter to bacteria (tiara_label)
    n0 = len(df)
    if args.filter_tiara_bacteria:
        df = df[df["tiara_label"].apply(tiara_is_bacteria)].copy()
        print(f"[INFO] Filter tiara_label == 'bacteria/archaea/prokarya': {n0} -> {len(df)} rows")

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

    # Map feature_id -> (contig, start, end) using contig+ordinal lookup in cds_by_contig
    mapped = 0
    unmapped = 0

    start_vals = []
    end_vals = []
    host_id_vals = []
    contig_vals = []

    for _, row in df.iterrows():
        src = (row.get("source", "") or "").strip().lower()
        fid = (row.get("feature_id", "") or "").strip()

        # default
        contig_vals.append(row.get("contig_id", "") or "")
        start_vals.append("")
        end_vals.append("")
        host_id_vals.append("")

        if not fid:
            unmapped += 1
            continue

        if src == "smorfinder":
            hit = smorf_index.get(fid)
            if not hit:
                unmapped += 1
                continue
            contig, s_i, e_i, strand = hit
            contig_vals[-1] = contig
            start_vals[-1] = str(s_i)
            end_vals[-1] = str(e_i)
            # host_id_vals stays blank here; hosts come from overlap step
            mapped += 1
            continue

        if src == "prodigal":

            # 1) First try direct lookup by feature_id (works for BYUMGR_* IDs)
            hit = cds_index.get(fid)
            if hit:
                contig_norm, s_i, e_i = hit
                contig_vals[-1] = contig_norm
                start_vals[-1] = str(s_i)
                end_vals[-1] = str(e_i)
                host_id_vals[-1] = fid
                mapped += 1
                continue

            # 2) Fallback: contig_###_ordinal style (if ever present)
            contig_key, ord_k = parse_feature_id(fid)
            if contig_key is None or ord_k is None:
                unmapped += 1
                continue

            contig_norm = normalize_contig_id(contig_key)
            cds_list = cds_by_contig.get(contig_norm)
            if not cds_list or ord_k < 1 or ord_k > len(cds_list):
                unmapped += 1
                continue

            cds = cds_list[ord_k - 1]
            contig_vals[-1] = contig_norm
            start_vals[-1] = str(cds["start"])
            end_vals[-1] = str(cds["end"])
            host_id_vals[-1] = str(cds.get("id", ""))
            mapped += 1
            continue

        # for "mixed" (or anything else), we intentionally do NOT assign coords
        unmapped += 1

    df["contig_id"] = contig_vals
    df["start"] = start_vals
    df["end"] = end_vals
    df["host_cds_id"] = host_id_vals

    print(f"[INFO] source-aware coord mapping: mapped={mapped}, unmapped={unmapped}")

    # ---------------------------------------------------------
    # SmORFinder-only overlap vs Prodigal CDS
    # Write results into EXISTING columns:
    #   - host_cds_id            (comma list of overlapping CDS IDs)
    #   - flag_overlap_fraction  (comma list of overlap fractions aligned to host_cds_id)
    #   - flag_embedded          (1 if ANY host fully contains smORF else 0)
    #   - host_cds_id_embedded   (comma list of host CDS IDs that fully contain smORF)
    # ---------------------------------------------------------
    ensure_cols(df, ["host_cds_id_embedded"])  # this is in your schema list
    # (host_cds_id, flag_overlap_fraction, flag_embedded already exist in your ensure_cols earlier)

    smorf_mask = df["source"].fillna("").astype(str).str.lower().eq("smorfinder")

    # Default-fill ONLY smorfinder rows; leave prodigal/mixed untouched
    df.loc[smorf_mask, "host_cds_id"] = ""
    df.loc[smorf_mask, "flag_overlap_fraction"] = ""
    df.loc[smorf_mask, "flag_embedded"] = ""
    df.loc[smorf_mask, "host_cds_id_embedded"] = ""

    # Iterate only smorfinder rows
    for idx, row in df.loc[smorf_mask].iterrows():
        contig = (row.get("contig_id", "") or "").strip()
        try:
            s_i = int(float((row.get("start", "") or "").strip()))
            e_i = int(float((row.get("end", "") or "").strip()))
        except Exception:
            # keep blanks if coords missing
            continue

        cds_list = cds_by_contig.get(contig, [])
        if not cds_list:
            # No Prodigal CDS on contig => embedded=0, overlap empty
            df.at[idx, "flag_embedded"] = "0"
            df.at[idx, "host_cds_id"] = ""
            df.at[idx, "flag_overlap_fraction"] = ""
            df.at[idx, "host_cds_id_embedded"] = ""
            continue

        host_ids, fracs, embedded_hosts = compute_overlaps_all_hosts(s_i, e_i, cds_list)

        # Write lists as comma-separated strings (aligned)
        # If you prefer semicolon, just change "," to ";"
        df.at[idx, "host_cds_id"] = ",".join([h for h in host_ids if h])
        df.at[idx, "flag_overlap_fraction"] = ",".join([f"{x:.3f}" for x in fracs])
        df.at[idx, "host_cds_id_embedded"] = ",".join([h for h in embedded_hosts if h])
        df.at[idx, "flag_embedded"] = "1" if len(embedded_hosts) > 0 else "0"

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

    # Step 2: embedded-in-longer-CDS (OPTIONAL, but DO NOT overwrite smorfinder)
    prodigal_mask = df["source"].fillna("").astype(str).str.lower().eq("prodigal")

    if args.prodigal_gff and len(cds_by_contig) > 0 and prodigal_mask.any():
        print("Applying Step 2 (prodigal-only): embedded-in-longer-CDS + overlap fraction...")

        step2 = compute_step2_embedded_and_overlap(
            df.loc[prodigal_mask].copy(),
            cds_by_contig=cds_by_contig,
            min_host_extra_bp=args.min_host_extra_bp
        )

        df.loc[prodigal_mask, "flag_embedded"] = step2["flag_embedded"]
        df.loc[prodigal_mask, "flag_overlap_fraction"] = step2["flag_overlap_fraction"]

        ensure_cols(df, ["host_cds_id_mapped", "host_cds_id_embedded"])
        df.loc[prodigal_mask, "host_cds_id_mapped"] = df.loc[prodigal_mask, "host_cds_id"]
        df.loc[prodigal_mask, "host_cds_id_embedded"] = step2["host_cds_id_step2"]
    else:
        print("[INFO] Step2 skipped (no prodigal rows or no CDS parsed).")

    df["confidence_tier"] = df["confidence_tier"].replace("", "UNASSESSED")

    print(f"[INFO] Rows: {len(df)}")
    print(f"[INFO] Missing contigs in bac_contigs.fasta (after normalization): {missing_contigs}")
    print(f"[INFO] Rows with non-integer coords: {bad_coords}")
    print(f"[INFO] flag_edge=1: {flagged}")

    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
