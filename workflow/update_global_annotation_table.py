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

# -----------------------------
# Annotation cutoffs
# -----------------------------

DIAMOND_MAX_EVALUE = 1e-5
DIAMOND_MIN_BITSCORE = 50.0
DIAMOND_MIN_PIDENT = 35.0
DIAMOND_MIN_QUERY_COV = 0.60

PFAM_MAX_IEVALUE = 1e-3
PFAM_MIN_DOMAIN_SCORE = 25.0
PFAM_MIN_DOMAIN_COV = 0.35

EGGNOG_MAX_EVALUE = 1e-3
EGGNOG_MIN_SCORE = 50.0


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

def safe_float(x, default=None):
    try:
        s = str(x).strip()
        if s in ("", "NA", "-", "nan", "None"):
            return default
        return float(s)
    except Exception:
        return default


def passes_diamond_filters(row) -> bool:
    evalue = safe_float(row.get("evalue"))
    bitscore = safe_float(row.get("bitscore"))
    pident = safe_float(row.get("pident"))
    aln_len = safe_float(row.get("aln_len"))
    qlen = safe_float(row.get("qlen"))

    if None in (evalue, bitscore, pident, aln_len, qlen):
        return False

    if qlen <= 0:
        return False

    query_cov = aln_len / qlen

    return (
        evalue <= DIAMOND_MAX_EVALUE
        and bitscore >= DIAMOND_MIN_BITSCORE
        and pident >= DIAMOND_MIN_PIDENT
        and query_cov >= DIAMOND_MIN_QUERY_COV
    )


def passes_eggnog_filters(row) -> bool:
    evalue = safe_float(row.get("evalue"))
    score = safe_float(row.get("score"))

    if evalue is None or score is None:
        return False

    return (
        evalue <= EGGNOG_MAX_EVALUE
        and score >= EGGNOG_MIN_SCORE
    )

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

    cols = [
        "qseqid", "sseqid", "stitle", "pident", "aln_len",
        "qlen", "slen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore"
    ]

    hits = pd.read_csv(args.diamond_tsv, sep="\t", names=cols, dtype=str, keep_default_na=False)

    if hits.empty:
        df = recompute_flags(df)
        write_table(df, args.table_tsv)
        return

    hits = hits[hits.apply(passes_diamond_filters, axis=1)].copy()

    if hits.empty:
        df = recompute_flags(df)
        write_table(df, args.table_tsv)
        return

    hits["query_cov"] = hits.apply(
        lambda r: safe_float(r["aln_len"], 0.0) / safe_float(r["qlen"], 1.0),
        axis=1,
    )

    hits["best_hit"] = hits.apply(
        lambda r: " | ".join([
            normalize_na(r["sseqid"]),
            normalize_na(r["stitle"]),
            f"pident={normalize_na(r['pident'])}",
            f"query_cov={r['query_cov']:.3f}",
            f"evalue={normalize_na(r['evalue'])}",
            f"bitscore={normalize_na(r['bitscore'])}",
        ]),
        axis=1,
    )

    hit_map = dict(zip(hits["qseqid"], hits["best_hit"]))

    df["best_ortholog_or_best_hit_description"] = df.apply(
        lambda r: hit_map.get(
            r["representative_id"],
            r["best_ortholog_or_best_hit_description"]
        ),
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
            hmm_len = safe_float(parts[2])
            query_id = parts[3]

            full_evalue = safe_float(parts[6])
            full_score = safe_float(parts[7])

            domain_i_eval = safe_float(parts[12])
            domain_score = safe_float(parts[13])

            hmm_from = safe_float(parts[15])
            hmm_to = safe_float(parts[16])

            if None in (hmm_len, domain_i_eval, domain_score, hmm_from, hmm_to):
                continue

            if hmm_len <= 0:
                continue

            domain_cov = (hmm_to - hmm_from + 1) / hmm_len

            passes_standard_domain_rule = (
                domain_i_eval <= PFAM_MAX_IEVALUE
                and domain_score >= PFAM_MIN_DOMAIN_SCORE
                and domain_cov >= PFAM_MIN_DOMAIN_COV
            )

            passes_strong_partial_rule = (
                domain_i_eval <= 1e-10
                and domain_score >= 40.0
            )

            if not (passes_standard_domain_rule or passes_strong_partial_rule):
                continue

            domain_match_type = (
                "standard"
                if passes_standard_domain_rule
                else "strong_partial"
            )

            label = (
                f"{domain_name}"
                f"({domain_acc};"
                f"i-evalue={domain_i_eval:.2e};"
                f"score={domain_score:.1f};"
                f"domain_cov={domain_cov:.3f};"
                f"match_type={domain_match_type})"
            )

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

    ann = ann[ann.apply(passes_eggnog_filters, axis=1)].copy()

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
        cols = ["qseqid", "sseqid", "stitle", "pident", "aln_len", "evalue", "bitscore"]
        hits = pd.read_csv(args.diamond_tsv, sep="\t", names=cols, dtype=str, keep_default_na=False)

        for _, row in hits.iterrows():
            qid = row["qseqid"]
            title = normalize_na(row["stitle"])

            # Try to extract organism from UniProt-style title, e.g. OS=Escherichia coli OX=...
            sp = "NA"
            m = re.search(r"OS=([^=]+?)(?:\sOX=|\sGN=|\sPE=|\sSV=|$)", title)
            if m:
                sp = m.group(1).strip()

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
