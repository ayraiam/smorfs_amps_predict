#!/usr/bin/env python3
import argparse
import math
import os
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

COLS = [
    "representative_id",
    "protein_length",
    "best_ortholog_or_best_hit_description",
    "conserved_domains",
    "enzyme_related_annotation",
    "go_terms",
    "kegg_ortholog_pathway_hint",
    "cog_functional_category",
    "taxonomic_hint",
    "annotation_confidence_tier",
    "hypothetical_or_unknown",
]

NA_COLS = [
    "best_ortholog_or_best_hit_description",
    "conserved_domains",
    "enzyme_related_annotation",
    "go_terms",
    "kegg_ortholog_pathway_hint",
    "cog_functional_category",
    "taxonomic_hint",
    "annotation_confidence_tier",
    "hypothetical_or_unknown",
]


def read_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in COLS:
        if c not in df.columns:
            df[c] = "NA"
    return df[COLS].copy()


def write_table(df: pd.DataFrame, path: str) -> None:
    df = df.copy()
    for c in COLS:
      if c not in df.columns:
        df[c] = "NA"
    df = df[COLS]
    df.to_csv(path, sep="\t", index=False)


def is_missing(x: str) -> bool:
    return x in ("", "NA", "-", None)


def normalize_na(x):
    if x is None:
        return "NA"
    if isinstance(x, float) and math.isnan(x):
        return "NA"
    s = str(x).strip()
    return "NA" if s == "" else s


def translate_record(seq: str) -> str:
    seq = re.sub(r"\s+", "", seq.upper())
    trim = len(seq) % 3
    if trim:
        seq = seq[: len(seq) - trim]
    aa = str(Seq(seq).translate(to_stop=False))
    aa = aa.rstrip("*")
    return aa


def recompute_flags(df: pd.DataFrame) -> pd.DataFrame:
    def _tier(row):
        orth = row["best_ortholog_or_best_hit_description"]
        dom = row["conserved_domains"]
        enz = row["enzyme_related_annotation"]
        go = row["go_terms"]
        keg = row["kegg_ortholog_pathway_hint"]
        cog = row["cog_functional_category"]
        tax = row["taxonomic_hint"]

        evidence_count = sum(not is_missing(v) for v in [orth, dom, enz, go, keg, cog, tax])

        text = " ".join([str(orth), str(enz), str(cog)]).lower()
        hypo = bool(re.search(r"\b(hypothetical|uncharacterized|unknown protein|predicted protein)\b", text))

        if evidence_count >= 4 and not hypo:
            return "high"
        if evidence_count >= 2 and not hypo:
            return "medium"
        if evidence_count >= 1:
            return "low"
        return "unknown"

    def _hypo(row):
        text = " ".join([
            str(row["best_ortholog_or_best_hit_description"]),
            str(row["enzyme_related_annotation"]),
            str(row["cog_functional_category"]),
        ]).lower()

        if re.search(r"\b(hypothetical|uncharacterized|unknown protein|predicted protein)\b", text):
            return "yes"

        if all(is_missing(row[c]) for c in [
            "best_ortholog_or_best_hit_description",
            "conserved_domains",
            "enzyme_related_annotation",
            "go_terms",
            "kegg_ortholog_pathway_hint",
            "cog_functional_category",
            "taxonomic_hint",
        ]):
            return "yes"
        return "no"

    df["annotation_confidence_tier"] = df.apply(_tier, axis=1)
    df["hypothetical_or_unknown"] = df.apply(_hypo, axis=1)
    return df


def cmd_init(args):
    records = []
    out_faa = Path(args.proteins_faa)
    out_faa.parent.mkdir(parents=True, exist_ok=True)

    with open(out_faa, "w") as fout:
        for rec in SeqIO.parse(args.input_fna, "fasta"):
            aa = translate_record(str(rec.seq))
            fout.write(f">{rec.id}\n{aa}\n")
            records.append({
                "representative_id": rec.id,
                "protein_length": str(len(aa)),
                **{c: "NA" for c in NA_COLS},
            })

    df = pd.DataFrame(records, columns=COLS)
    df = recompute_flags(df)
    write_table(df, args.output_tsv)


def cmd_merge_orthology(args):
    df = read_table(args.table_tsv)
    if not Path(args.diamond_tsv).exists():
        raise FileNotFoundError(args.diamond_tsv)

    cols = ["qseqid", "sseqid", "stitle", "pident", "aln_len", "evalue", "bitscore", "sscinames"]
    hits = pd.read_csv(args.diamond_tsv, sep="\t", names=cols, dtype=str, keep_default_na=False)
    if hits.empty:
        df = recompute_flags(df)
        write_table(df, args.table_tsv)
        return

    hits["best_hit"] = hits.apply(
        lambda r: " | ".join([
            normalize_na(r["sseqid"]),
            normalize_na(r["stitle"]),
            f"pident={normalize_na(r['pident'])}",
            f"evalue={normalize_na(r['evalue'])}",
            f"bitscore={normalize_na(r['bitscore'])}",
        ]),
        axis=1,
    )
    hit_map = dict(zip(hits["qseqid"], hits["best_hit"]))

    df["best_ortholog_or_best_hit_description"] = df.apply(
        lambda r: hit_map.get(r["representative_id"], r["best_ortholog_or_best_hit_description"]),
        axis=1,
    )

    df = recompute_flags(df)
    write_table(df, args.table_tsv)


def parse_domtblout(path: str):
    parsed = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 23:
                continue
            domain_name = parts[0]
            domain_acc = parts[1]
            query_id = parts[3]
            full_evalue = parts[6]
            dom_score = parts[13]
            label = f"{domain_name}({domain_acc};e={full_evalue};score={dom_score})"
            parsed.setdefault(query_id, [])
            if label not in parsed[query_id]:
                parsed[query_id].append(label)
    return {k: "; ".join(v[:5]) for k, v in parsed.items()}


def cmd_merge_domains(args):
    df = read_table(args.table_tsv)
    if not Path(args.domtblout).exists():
        raise FileNotFoundError(args.domtblout)

    domain_map = parse_domtblout(args.domtblout)
    df["conserved_domains"] = df.apply(
        lambda r: domain_map.get(r["representative_id"], r["conserved_domains"]),
        axis=1,
    )
    df = recompute_flags(df)
    write_table(df, args.table_tsv)


def load_eggnog(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            rows.append(line.rstrip("\n").split("\t"))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows, columns=header)
    return df


def first_existing(df: pd.DataFrame, names):
    for n in names:
        if n in df.columns:
            return n
    return None


def cmd_merge_functional(args):
    df = read_table(args.table_tsv)
    if not Path(args.eggnog_annotations).exists():
        raise FileNotFoundError(args.eggnog_annotations)

    ann = load_eggnog(args.eggnog_annotations)
    if ann.empty:
        df = recompute_flags(df)
        write_table(df, args.table_tsv)
        return

    qcol = first_existing(ann, ["query", "#query"])
    desc_col = first_existing(ann, ["Description", "description"])
    ec_col = first_existing(ann, ["EC"])
    go_col = first_existing(ann, ["GOs"])
    kegg_ko_col = first_existing(ann, ["KEGG_ko"])
    kegg_pw_col = first_existing(ann, ["KEGG_Pathway"])
    cog_col = first_existing(ann, ["COG_category"])

    merged = {}
    for _, row in ann.iterrows():
        qid = row[qcol]
        desc = normalize_na(row[desc_col]) if desc_col else "NA"
        ec = normalize_na(row[ec_col]) if ec_col else "NA"
        gos = normalize_na(row[go_col]) if go_col else "NA"
        kos = normalize_na(row[kegg_ko_col]) if kegg_ko_col else "NA"
        pws = normalize_na(row[kegg_pw_col]) if kegg_pw_col else "NA"
        cog = normalize_na(row[cog_col]) if cog_col else "NA"

        kegg_hint = "NA"
        parts = []
        if not is_missing(kos):
            parts.append(f"KO={kos}")
        if not is_missing(pws):
            parts.append(f"Pathway={pws}")
        if parts:
            kegg_hint = " | ".join(parts)

        merged[qid] = {
            "desc": desc,
            "ec": ec,
            "gos": gos,
            "kegg_hint": kegg_hint,
            "cog": cog,
        }

    def update_row(row):
        info = merged.get(row["representative_id"])
        if info is None:
            return row

        if is_missing(row["best_ortholog_or_best_hit_description"]) and not is_missing(info["desc"]):
            row["best_ortholog_or_best_hit_description"] = info["desc"]
        if not is_missing(info["ec"]):
            row["enzyme_related_annotation"] = info["ec"]
        if not is_missing(info["gos"]):
            row["go_terms"] = info["gos"]
        if not is_missing(info["kegg_hint"]):
            row["kegg_ortholog_pathway_hint"] = info["kegg_hint"]
        if not is_missing(info["cog"]):
            row["cog_functional_category"] = info["cog"]
        return row

    df = df.apply(update_row, axis=1)
    df = recompute_flags(df)
    write_table(df, args.table_tsv)


def cmd_merge_taxonomy(args):
    df = read_table(args.table_tsv)

    tax_map = {}

    if args.diamond_tsv and Path(args.diamond_tsv).exists():
        cols = ["qseqid", "sseqid", "stitle", "pident", "aln_len", "evalue", "bitscore", "sscinames"]
        hits = pd.read_csv(args.diamond_tsv, sep="\t", names=cols, dtype=str, keep_default_na=False)
        for _, row in hits.iterrows():
            qid = row["qseqid"]
            sp = normalize_na(row["sscinames"])
            if not is_missing(sp):
                tax_map[qid] = sp

    if args.eggnog_annotations and Path(args.eggnog_annotations).exists():
        ann = load_eggnog(args.eggnog_annotations)
        if not ann.empty:
            qcol = first_existing(ann, ["query", "#query"])
            tax_scope_col = first_existing(ann, ["tax_scope"])
            max_annot_col = first_existing(ann, ["max_annot_lvl"])
            for _, row in ann.iterrows():
                qid = row[qcol]
                if qid in tax_map:
                    continue
                parts = []
                if tax_scope_col and not is_missing(normalize_na(row[tax_scope_col])):
                    parts.append(f"tax_scope={normalize_na(row[tax_scope_col])}")
                if max_annot_col and not is_missing(normalize_na(row[max_annot_col])):
                    parts.append(f"max_annot_lvl={normalize_na(row[max_annot_col])}")
                if parts:
                    tax_map[qid] = " | ".join(parts)

    df["taxonomic_hint"] = df.apply(
        lambda r: tax_map.get(r["representative_id"], r["taxonomic_hint"]),
        axis=1,
    )

    df = recompute_flags(df)
    write_table(df, args.table_tsv)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("init")
    p.add_argument("--input-fna", required=True)
    p.add_argument("--proteins-faa", required=True)
    p.add_argument("--output-tsv", required=True)
    p.set_defaults(func=cmd_init)

    p = sub.add_parser("merge-orthology")
    p.add_argument("--table-tsv", required=True)
    p.add_argument("--diamond-tsv", required=True)
    p.set_defaults(func=cmd_merge_orthology)

    p = sub.add_parser("merge-domains")
    p.add_argument("--table-tsv", required=True)
    p.add_argument("--domtblout", required=True)
    p.set_defaults(func=cmd_merge_domains)

    p = sub.add_parser("merge-functional")
    p.add_argument("--table-tsv", required=True)
    p.add_argument("--eggnog-annotations", required=True)
    p.set_defaults(func=cmd_merge_functional)

    p = sub.add_parser("merge-taxonomy")
    p.add_argument("--table-tsv", required=True)
    p.add_argument("--diamond-tsv", required=False, default="")
    p.add_argument("--eggnog-annotations", required=False, default="")
    p.set_defaults(func=cmd_merge_taxonomy)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
