#!/usr/bin/env python3
import argparse
from collections import defaultdict
import pandas as pd
from Bio import SeqIO

def ensure_cols(df, cols):
    for c in cols:
        if c not in df.columns:
            df[c] = ""

def build_contig_len_dict(fasta_path):
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}

def parse_metaeuk_gff(gff_path):
    feats_by_contig = defaultdict(list)
    feat_index = {}

    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = fields

            if feature_type != "CDS":
                continue

            ad = {}
            for attr in attrs.split(";"):
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    ad[k.strip()] = v.strip()

            fid = ad.get("TCS_ID", "") or ad.get("ID", "")
            if not fid:
                continue

            start_i = int(start)
            end_i = int(end)

            feats_by_contig[seqid].append({
                "id": fid,
                "start": start_i,
                "end": end_i,
                "strand": strand,
                "feature_type": feature_type
            })
            feat_index[fid] = (seqid, start_i, end_i, strand)

    return feats_by_contig, feat_index

def parse_metaeuk_feature_id(fid: str):
    """
    MetaEuk feature_id format in predicted_smorfs.tsv:
      Target|contig|strand|score|evalue|exon_count|start|end|exon_blocks

    Example:
      UniRef50_A0A820GSM0|contig_6867|-|74|1.83e-12|1|4211|4321|4321[4321]:4211[4211]:111[111]

    Returns:
      contig_id, strand, start, end
    """
    fid = (fid or "").strip()
    parts = fid.split("|")
    if len(parts) < 8:
        return None, None, None, None

    try:
        contig = parts[1]
        strand = parts[2]
        start = int(parts[6])
        end = int(parts[7])
        return contig, strand, start, end
    except Exception:
        return None, None, None, None

def load_mmseqs_cluster_map(cluster_map_tsv: str) -> pd.DataFrame:
    cm = pd.read_csv(cluster_map_tsv, sep="\t", header=None, dtype=str).fillna("")
    if cm.shape[1] < 2:
        raise SystemExit(f"ERROR: cluster_map.tsv must have at least 2 columns. Got {cm.shape[1]}")
    cm = cm.iloc[:, :2].copy()
    cm.columns = ["member_id", "rep_id"]
    return cm

def attach_cluster_stats(df: pd.DataFrame, cluster_map_tsv: str, env_label: str | None = None) -> pd.DataFrame:
    cm = load_mmseqs_cluster_map(cluster_map_tsv)

    df = df.copy()
    if env_label:
        df["member_id"] = env_label + "|" + df["feature_id"].astype(str)
    else:
        df["member_id"] = df["feature_id"].astype(str)

    to_drop = ["cluster_id", "cluster_size", "cluster_env_count",
               "flag_cluster_recurrent", "flag_cross_environment"]
    df.drop(columns=[c for c in to_drop if c in df.columns], inplace=True, errors="ignore")

    rep_size = cm.groupby("rep_id")["member_id"].size().rename("cluster_size").reset_index()
    cm["env_global"] = cm["member_id"].str.split("|", n=1).str[0] + "_GLOBAL"
    rep_env = cm.groupby("rep_id")["env_global"].nunique().rename("cluster_env_count").reset_index()
    rep_stats = rep_size.merge(rep_env, on="rep_id", how="left")

    member_to_rep = cm[["member_id", "rep_id"]].drop_duplicates()

    out = df.merge(member_to_rep, on="member_id", how="left")
    out = out.merge(rep_stats, on="rep_id", how="left")

    out["cluster_id"] = out["rep_id"].fillna("")
    out["cluster_size"] = out["cluster_size"].fillna("").astype(str)
    out["cluster_env_count"] = out["cluster_env_count"].fillna("").astype(str)

    cs = pd.to_numeric(out["cluster_size"], errors="coerce").fillna(0).astype(int)
    ec = pd.to_numeric(out["cluster_env_count"], errors="coerce").fillna(0).astype(int)

    out["flag_cluster_recurrent"] = (cs >= 2).astype(int).astype(str)
    out["flag_cross_environment"] = (ec >= 2).astype(int).astype(str)

    out.drop(columns=["rep_id"], inplace=True, errors="ignore")
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-tsv", required=True)
    ap.add_argument("--metaeuk-gff", required=True)
    ap.add_argument("--fungi-contigs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--edge-bp", type=int, default=50)
    ap.add_argument("--run-step1", type=int, choices=[0,1], default=1)
    ap.add_argument("--run-step2", type=int, choices=[0,1], default=1)
    ap.add_argument("--run-step3", type=int, choices=[0,1], default=1)
    ap.add_argument("--cluster-map", default="")
    ap.add_argument("--env-label", default="")
    ap.add_argument(
        "--cluster-only",
        action="store_true",
        help="Only attach MMseqs cluster stats (Step 3) and write output; skip Steps 1/2.",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str).copy()
    ensure_cols(df, [
        "start", "end", "flag_edge", "dist_left", "dist_right",
        "flag_embedded", "flag_overlap_fraction", "host_cds_id",
        "host_cds_id_embedded", "confidence_tier",
        "cluster_id", "cluster_size", "cluster_env_count",
        "flag_cluster_recurrent", "flag_cross_environment"
    ])

    contig_len = build_contig_len_dict(args.fungi_contigs)
    metaeuk_by_contig, _ = parse_metaeuk_gff(args.metaeuk_gff)

    keep = df["tiara_label"].fillna("").str.lower().str.contains("euk|fung")
    df = df.loc[keep].copy()

    if args.cluster_only:
        if args.cluster_map:
            env_label = args.env_label.strip() or None
            df = attach_cluster_stats(df, args.cluster_map, env_label=env_label)
            print("[INFO] (cluster-only) Step3 cluster stats attached.")
        else:
            print("[INFO] (cluster-only) requested but no cluster_map provided; skipped.")
        df["confidence_tier"] = df["confidence_tier"].replace("", "UNASSESSED")
        df.to_csv(args.out, sep="\t", index=False)
        print(f"[OK] (cluster-only) wrote: {args.out}")
        raise SystemExit(0)

    # STEP 1
    if args.run_step1 == 1:
        mapped = 0
        unmapped = 0

        for idx, row in df.iterrows():
            src = str(row.get("source", "")).strip().lower()
            fid = str(row.get("feature_id", "")).strip()

            if src != "metaeuk":
                unmapped += 1
                continue

            contig, strand, s, e = parse_metaeuk_feature_id(fid)
            if contig is None:
                unmapped += 1
                continue

            df.at[idx, "contig_id"] = contig
            df.at[idx, "start"] = str(s)
            df.at[idx, "end"] = str(e)
            mapped += 1

        for idx, row in df.iterrows():
            contig = str(row.get("contig_id", "")).strip()
            try:
                s = int(float(str(row.get("start", "")).strip()))
                e = int(float(str(row.get("end", "")).strip()))
            except Exception:
                continue

            L = contig_len.get(contig)
            if L is None:
                continue

            dl = s - 1
            dr = L - e
            df.at[idx, "dist_left"] = str(dl)
            df.at[idx, "dist_right"] = str(dr)
            df.at[idx, "flag_edge"] = "1" if (dl < args.edge_bp or dr < args.edge_bp) else "0"

        print(f"[INFO] Step1 mapped={mapped}, unmapped={unmapped}")
    else:
        print("[INFO] Step1 skipped.")

    # STEP 2
    if args.run_step2 == 1:
        print("[INFO] Applying Step2: MetaEuk overlap / embedded checks")
        for idx, row in df.iterrows():
            contig = str(row.get("contig_id", "")).strip()
            try:
                s = int(float(str(row.get("start", "")).strip()))
                e = int(float(str(row.get("end", "")).strip()))
            except Exception:
                continue

            cds_list = metaeuk_by_contig.get(contig, [])
            if not cds_list:
                df.at[idx, "flag_embedded"] = "0"
                df.at[idx, "flag_overlap_fraction"] = ""
                df.at[idx, "host_cds_id"] = ""
                df.at[idx, "host_cds_id_embedded"] = ""
                continue

            sm_len = e - s + 1
            host_ids = []
            fracs = []
            embedded_hosts = []

            for cds in cds_list:
                cs = int(cds["start"])
                ce = int(cds["end"])
                ov = max(0, min(e, ce) - max(s, cs) + 1)
                if ov <= 0:
                    continue

                frac = ov / sm_len if sm_len > 0 else 0.0
                host_ids.append(str(cds.get("id", "")))
                fracs.append(frac)

                if cs <= s and ce >= e:
                    embedded_hosts.append(str(cds.get("id", "")))

            df.at[idx, "host_cds_id"] = ",".join([x for x in host_ids if x])
            df.at[idx, "flag_overlap_fraction"] = ",".join([f"{x:.3f}" for x in fracs])
            df.at[idx, "host_cds_id_embedded"] = ",".join([x for x in embedded_hosts if x])
            df.at[idx, "flag_embedded"] = "1" if embedded_hosts else "0"
    else:
        print("[INFO] Step2 skipped.")

    # STEP 3
    if args.run_step3 == 1:
        if args.cluster_map:
            env_label = args.env_label.strip() or None
            df = attach_cluster_stats(df, args.cluster_map, env_label=env_label)
            print("[INFO] Step3 cluster stats attached.")
        else:
            print("[INFO] Step3 requested but no cluster_map provided; skipped.")
    else:
        print("[INFO] Step3 skipped.")

    df["confidence_tier"] = df["confidence_tier"].replace("", "UNASSESSED")
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote: {args.out}")

if __name__ == "__main__":
    main()
