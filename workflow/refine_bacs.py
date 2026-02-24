#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict


def build_contig_len_dict(fasta_path):
    contig_len = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        contig_len[record.id] = len(record.seq)
    return contig_len


def build_cds_dict_from_gff(gff_path):
    cds_by_contig = defaultdict(list)

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields

            if feature_type != "CDS":
                continue

            start = int(start)
            end = int(end)

            # Extract ID if present
            cds_id = ""
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    cds_id = attr.replace("ID=", "")
                    break

            cds_by_contig[seqid].append(
                {
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "id": cds_id
                }
            )

    return cds_by_contig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-tsv", required=True)
    ap.add_argument("--prodigal-gff", required=True)
    ap.add_argument("--bac-contigs", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpus", type=int, default=8)
    args = ap.parse_args()

    print("=======================================")
    print(f"Sample: {args.sample}")
    print("Building initial dictionaries...")
    print("=======================================")

    # Contig lengths
    contig_len = build_contig_len_dict(args.bac_contigs)
    print(f"[INFO] Contigs loaded: {len(contig_len)}")

    # Print first 5
    for i, (k, v) in enumerate(contig_len.items()):
        print(f"  {k} -> length={v}")
        if i == 4:
            break

    # CDS intervals
    cds_by_contig = build_cds_dict_from_gff(args.prodigal_gff)
    print(f"[INFO] Contigs with CDS: {len(cds_by_contig)}")

    # Print first 3 contigs with CDS
    for i, (contig, cds_list) in enumerate(cds_by_contig.items()):
        print(f"\n  {contig}: {len(cds_list)} CDS entries")
        for j, cds in enumerate(cds_list[:3]):
            print(f"    {cds}")
        if i == 2:
            break

    print("\n[OK] Initial dicts built successfully.")
    print("Refinement logic not yet applied.")

    # Just copy input to output for now
    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote: {args.out}")


if __name__ == "__main__":
    main()
