#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description="Build GLOBAL representative CDS FASTA from MMseqs cluster_map.tsv"
    )
    ap.add_argument("--cluster-map", required=True,
                    help="Path to cluster_map.tsv (member<TAB>rep)")
    ap.add_argument("--results-dir", required=True,
                    help="Base results dir, e.g. results")
    ap.add_argument("--out-fasta", required=True,
                    help="Output FASTA of representative CDSs")
    ap.add_argument("--out-meta", required=True,
                    help="Output metadata TSV for representative CDSs")
    return ap.parse_args()


def read_unique_reps(cluster_map_path):
    reps = set()
    with open(cluster_map_path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if len(row) < 2:
                raise ValueError(
                    f"cluster_map.tsv must have at least 2 columns; got: {row}"
                )
            rep = row[1].strip()
            if rep:
                reps.add(rep)
    return sorted(reps)


def split_rep_id(rep_id):
    if "|" not in rep_id:
        raise ValueError(f"Representative ID lacks '|': {rep_id}")
    env, cds_id = rep_id.split("|", 1)
    if not env or not cds_id:
        raise ValueError(f"Malformed representative ID: {rep_id}")
    return env, cds_id


def fasta_iter(path):
    header = None
    seq_chunks = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    args = parse_args()

    cluster_map = Path(args.cluster_map)
    results_dir = Path(args.results_dir)
    out_fasta = Path(args.out_fasta)
    out_meta = Path(args.out_meta)

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_meta.parent.mkdir(parents=True, exist_ok=True)

    reps = read_unique_reps(cluster_map)

    # group wanted CDS IDs by environment
    wanted_by_env = defaultdict(set)
    for rep_id in reps:
        env, cds_id = split_rep_id(rep_id)
        wanted_by_env[env].add(cds_id)

    found = set()
    written = 0

    with open(out_fasta, "w") as fasta_out, open(out_meta, "w", newline="") as meta_out:
        meta_writer = csv.writer(meta_out, delimiter="\t")
        meta_writer.writerow([
            "rep_id",
            "environment",
            "cds_id",
            "source_fasta",
            "original_header"
        ])

        for env, wanted_ids in sorted(wanted_by_env.items()):
            src_fasta = results_dir / "smorfs" / f"{env}_GLOBAL" / "catalog" / "cds_all.fna"
            if not src_fasta.is_file():
                raise FileNotFoundError(f"Missing source FASTA: {src_fasta}")

            for header, seq in fasta_iter(src_fasta):
                first_token = header.split(" ", 1)[0]
                if first_token in wanted_ids:
                    rep_id = f"{env}|{first_token}"
                    rest = header[len(first_token):]  # preserve metadata after first token
                    fasta_out.write(f">{rep_id}{rest}\n")
                    for i in range(0, len(seq), 80):
                        fasta_out.write(seq[i:i+80] + "\n")

                    meta_writer.writerow([
                        rep_id,
                        env,
                        first_token,
                        str(src_fasta),
                        header
                    ])
                    found.add(rep_id)
                    written += 1

    missing = [r for r in reps if r not in found]
    if missing:
        preview = "\n".join(missing[:20])
        raise RuntimeError(
            f"Could not find {len(missing)} representative CDS IDs in source cds_all.fna files.\n"
            f"First missing IDs:\n{preview}"
        )

    print(f"[INFO] Wrote {written} representative CDS sequences to: {out_fasta}", file=sys.stderr)
    print(f"[INFO] Metadata TSV: {out_meta}", file=sys.stderr)


if __name__ == "__main__":
    main()
