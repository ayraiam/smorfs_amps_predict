#!/usr/bin/env python3

import argparse
import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx


def parse_args():
    ap = argparse.ArgumentParser(
        description="Build global cross-layer annotation network: Pfam, KO, GO, COG, EC"
    )
    ap.add_argument("--annotation-file", required=True)
    ap.add_argument("--count-matrix", required=True)
    ap.add_argument("--metadata-file", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--feature-col", default="representative_id")
    ap.add_argument("--min-edge-support", type=int, default=2)
    ap.add_argument("--min-edge-sum-cpm", type=float, default=0.0)
    ap.add_argument("--top-edge-percentile", type=float, default=0.0,
                    help="Optional percentile cutoff for sum_global_mean_cpm, e.g. 75")
    ap.add_argument(
    "--step",
    default="all",
    choices=[
        "all",
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "cds-universe",
        "clean-annotations",
        "abundance",
        "merge",
        "raw-edges",
        "collapse-edges",
        "nodes",
        "filter",
        "communities",
        "export",
    ],
    help="Network step to run"
    )
    return ap.parse_args()


def norm_na(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if x.upper() in {"NA", "N/A", "NONE", "NULL", ".", ""}:
        return ""
    return x


def split_semicolon_comma(x):
    x = norm_na(x)
    if not x:
        return []
    parts = re.split(r"[;,]\s*", x)
    return [p.strip() for p in parts if p.strip()]


def parse_pfams(x):
    vals = []
    for item in split_semicolon_comma(x):
        m = re.search(r"(PF\d{5})", item)
        if m:
            pfid = m.group(1)
            label = item.split("(")[0].strip() if "(" in item else item.strip()
            vals.append(("Pfam", pfid, label, f"PFAM:{pfid}"))
    return vals


def parse_kos(x):
    vals = []
    x = norm_na(x)
    if not x:
        return vals

    for ko in sorted(set(re.findall(r"K\d{5}", x))):
        vals.append(("KO", ko, ko, f"KO:{ko}"))
    return vals


def parse_gos(x):
    vals = []
    x = norm_na(x)
    if not x:
        return vals

    for go in sorted(set(re.findall(r"GO:\d{7}", x))):
        vals.append(("GO", go, go, f"GO:{go}"))
    return vals


def parse_ecs(x):
    vals = []
    x = norm_na(x)
    if not x:
        return vals

    # captures EC:1.2.3.4 or 1.2.3.4
    ecs = re.findall(r"(?:EC[:=]?)?\b(\d+\.\d+\.\d+\.\d+)\b", x)
    for ec in sorted(set(ecs)):
        vals.append(("EC", ec, ec, f"EC:{ec}"))
    return vals


def parse_cogs(x):
    vals = []
    x = norm_na(x)
    if not x:
        return vals

    # COG categories in your table are often letters, sometimes multiple letters
    candidates = re.findall(r"\b[A-Z]\b", x)
    if not candidates and len(x) <= 5:
        candidates = list(x)

    for cog in sorted(set(candidates)):
        if re.match(r"^[A-Z]$", cog):
            vals.append(("COG", cog, cog, f"COG:{cog}"))
    return vals


def find_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def compute_cds_abundance(count_file):
    counts = pd.read_csv(count_file, sep="\t")

    feature_col = counts.columns[0]
    sample_cols = counts.columns[1:]

    numeric = counts[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    lib_sizes = numeric.sum(axis=0)
    lib_sizes = lib_sizes.replace(0, np.nan)

    cpm = numeric.div(lib_sizes, axis=1) * 1_000_000

    out = pd.DataFrame({
        "representative_id": counts[feature_col].astype(str),
        "total_count": numeric.sum(axis=1),
        "global_mean_cpm": cpm.mean(axis=1),
        "global_median_cpm": cpm.median(axis=1),
        "global_prevalence": (numeric > 0).sum(axis=1) / numeric.shape[1],
        "n_samples_detected": (numeric > 0).sum(axis=1),
        "n_samples_total": numeric.shape[1],
    })

    return out


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def should_run(step_name):
        aliases = {
            "1": "cds-universe",
            "2": "clean-annotations",
            "3": "abundance",
            "4": "merge",
            "5": "raw-edges",
            "6": "collapse-edges",
            "7": "nodes",
            "8": "filter",
            "9": "communities",
            "10": "export",
        }

        if args.step == "all":
            return True

        wanted = aliases.get(args.step, args.step)
        return wanted == step_name

    def require_file(path):
        path = Path(path)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Required input file missing or empty: {path}")
        return path

    # -------------------------
    # Step 1 — Define CDS universe
    # -------------------------
    if should_run("cds-universe"):
        print("[INFO] Running Step 1 — Define CDS universe")

        annot = pd.read_csv(args.annotation_file, sep="\t")

        if args.feature_col not in annot.columns:
            raise ValueError(
                f"Feature column '{args.feature_col}' not found. "
                f"Available columns: {list(annot.columns)}"
            )

        base_cols = [args.feature_col]
        for c in [
            "protein_length",
            "annotation_confidence_tier",
            "hypothetical_or_unknown",
            "taxonomic_hint",
        ]:
            if c in annot.columns:
                base_cols.append(c)

        cds_universe = annot[base_cols].drop_duplicates()
        cds_universe = cds_universe.rename(columns={args.feature_col: "representative_id"})
        cds_universe.to_csv(outdir / "01_cds_universe.tsv", sep="\t", index=False)

    # -------------------------
    # Step 2 — Clean annotation columns
    # -------------------------
    if should_run("clean-annotations"):
        print("[INFO] Running Step 2 — Clean annotation columns")

        annot = pd.read_csv(args.annotation_file, sep="\t")

        if args.feature_col not in annot.columns:
            raise ValueError(
                f"Feature column '{args.feature_col}' not found. "
                f"Available columns: {list(annot.columns)}"
            )

        pfam_col = find_col(annot, ["conserved_domains", "pfam", "pfam_domains"])
        ko_col = find_col(annot, ["kegg_ortholog_pathway_hint", "kegg", "ko", "kegg_ko"])
        go_col = find_col(annot, ["go_terms", "go", "gos"])
        cog_col = find_col(annot, ["cog_functional_category", "cog", "cog_category"])
        ec_col = find_col(annot, ["enzyme_related_annotation", "ec", "ecs", "enzyme"])

        long_rows = []

        for _, row in annot.iterrows():
            rid = str(row[args.feature_col])

            anns = []
            if pfam_col:
                anns.extend(parse_pfams(row.get(pfam_col, "")))
            if ko_col:
                anns.extend(parse_kos(row.get(ko_col, "")))
            if go_col:
                anns.extend(parse_gos(row.get(go_col, "")))
            if cog_col:
                anns.extend(parse_cogs(row.get(cog_col, "")))
            if ec_col:
                anns.extend(parse_ecs(row.get(ec_col, "")))

            seen = set()
            for annotation_type, annotation_id, annotation_label, node_id in anns:
                key = (rid, node_id)
                if key in seen:
                    continue
                seen.add(key)

                long_rows.append({
                    "representative_id": rid,
                    "annotation_type": annotation_type,
                    "annotation_id": annotation_id,
                    "annotation_label": annotation_label,
                    "node_id": node_id,
                })

        annotation_long = pd.DataFrame(long_rows)
        annotation_long.to_csv(outdir / "02_clean_annotation_long.tsv", sep="\t", index=False)

    # -------------------------
    # Step 3 — Add abundance information
    # -------------------------
    if should_run("abundance"):
        print("[INFO] Running Step 3 — Add abundance information")

        abundance = compute_cds_abundance(args.count_matrix)
        abundance.to_csv(outdir / "03_cds_abundance_summary.tsv", sep="\t", index=False)

    # -------------------------
    # Step 4 — Merge annotation + abundance
    # -------------------------
    if should_run("merge"):
        print("[INFO] Running Step 4 — Merge annotation + abundance")

        annotation_long = pd.read_csv(
            require_file(outdir / "02_clean_annotation_long.tsv"),
            sep="\t"
        )
        abundance = pd.read_csv(
            require_file(outdir / "03_cds_abundance_summary.tsv"),
            sep="\t"
        )

        long_abund = annotation_long.merge(abundance, on="representative_id", how="left")
        long_abund.to_csv(outdir / "04_annotation_long_with_abundance.tsv", sep="\t", index=False)

    # -------------------------
    # Step 5 — Create pairwise annotation edges per CDS
    # -------------------------
    if should_run("raw-edges"):
        print("[INFO] Running Step 5 — Create pairwise annotation edges per CDS")

        long_abund = pd.read_csv(
            require_file(outdir / "04_annotation_long_with_abundance.tsv"),
            sep="\t"
        )
        abundance = pd.read_csv(
            require_file(outdir / "03_cds_abundance_summary.tsv"),
            sep="\t"
        )

        edge_rows = []

        for rid, sub in long_abund.groupby("representative_id"):
            nodes = sub[["node_id", "annotation_type", "annotation_label"]].drop_duplicates()

            if nodes.shape[0] < 2:
                continue

            abundance_row = abundance[abundance["representative_id"] == rid]
            if abundance_row.empty:
                cds_mean_cpm = 0.0
                cds_median_cpm = 0.0
                cds_prev = 0.0
            else:
                cds_mean_cpm = float(abundance_row["global_mean_cpm"].iloc[0])
                cds_median_cpm = float(abundance_row["global_median_cpm"].iloc[0])
                cds_prev = float(abundance_row["global_prevalence"].iloc[0])

            records = nodes.to_dict("records")

            for a, b in itertools.combinations(records, 2):
                s = a["node_id"]
                t = b["node_id"]

                if s == t:
                    continue

                if s > t:
                    a, b = b, a
                    s, t = t, s

                edge_layer = f"{a['annotation_type']}-{b['annotation_type']}"

                edge_rows.append({
                    "representative_id": rid,
                    "source_node": s,
                    "target_node": t,
                    "source_type": a["annotation_type"],
                    "target_type": b["annotation_type"],
                    "source_label": a["annotation_label"],
                    "target_label": b["annotation_label"],
                    "edge_layer": edge_layer,
                    "cds_global_mean_cpm": cds_mean_cpm,
                    "cds_global_median_cpm": cds_median_cpm,
                    "cds_global_prevalence": cds_prev,
                })

        raw_edges = pd.DataFrame(edge_rows)
        raw_edges.to_csv(outdir / "05_raw_edges_by_cds.tsv", sep="\t", index=False)

    # -------------------------
    # Step 6 — Collapse repeated edges
    # -------------------------
    if should_run("collapse-edges"):
        print("[INFO] Running Step 6 — Collapse repeated edges")

        raw_edges = pd.read_csv(
            require_file(outdir / "05_raw_edges_by_cds.tsv"),
            sep="\t"
        )

        if raw_edges.empty:
            raise RuntimeError("No edges were created. Check annotation parsing.")

        grouped = (
            raw_edges
            .groupby(
                [
                    "source_node",
                    "target_node",
                    "source_type",
                    "target_type",
                    "source_label",
                    "target_label",
                    "edge_layer",
                ],
                as_index=False,
            )
            .agg(
                n_cds_supporting_edge=("representative_id", "nunique"),
                sum_global_mean_cpm=("cds_global_mean_cpm", "sum"),
                median_global_mean_cpm=("cds_global_mean_cpm", "median"),
                mean_global_prevalence=("cds_global_prevalence", "mean"),
                supporting_cds_list=(
                    "representative_id",
                    lambda x: ";".join(sorted(set(map(str, x))))
                ),
            )
        )

        grouped.to_csv(outdir / "06_global_network_edges.tsv", sep="\t", index=False)

    # -------------------------
    # Step 7 — Create node table
    # -------------------------
    if should_run("nodes"):
        print("[INFO] Running Step 7 — Create node table")

        long_abund = pd.read_csv(
            require_file(outdir / "04_annotation_long_with_abundance.tsv"),
            sep="\t"
        )
        grouped = pd.read_csv(
            require_file(outdir / "06_global_network_edges.tsv"),
            sep="\t"
        )

        node_base = (
            long_abund
            .groupby(
                ["node_id", "annotation_type", "annotation_id", "annotation_label"],
                as_index=False,
            )
            .agg(
                n_cds_with_node=("representative_id", "nunique"),
                sum_global_mean_cpm=("global_mean_cpm", "sum"),
                median_global_mean_cpm=("global_mean_cpm", "median"),
                mean_global_prevalence=("global_prevalence", "mean"),
            )
        )

        G_full = nx.Graph()

        for _, r in node_base.iterrows():
            G_full.add_node(
                r["node_id"],
                annotation_type=r["annotation_type"],
                annotation_id=r["annotation_id"],
                annotation_label=r["annotation_label"],
                n_cds_with_node=int(r["n_cds_with_node"]),
                sum_global_mean_cpm=float(r["sum_global_mean_cpm"]),
                median_global_mean_cpm=float(r["median_global_mean_cpm"]),
                mean_global_prevalence=float(r["mean_global_prevalence"]),
            )

        for _, r in grouped.iterrows():
            G_full.add_edge(
                r["source_node"],
                r["target_node"],
                edge_layer=r["edge_layer"],
                n_cds_supporting_edge=int(r["n_cds_supporting_edge"]),
                sum_global_mean_cpm=float(r["sum_global_mean_cpm"]),
                median_global_mean_cpm=float(r["median_global_mean_cpm"]),
                mean_global_prevalence=float(r["mean_global_prevalence"]),
            )

        degree = dict(G_full.degree())
        node_base["degree"] = node_base["node_id"].map(degree).fillna(0).astype(int)
        node_base.to_csv(outdir / "07_global_network_nodes.tsv", sep="\t", index=False)

    # -------------------------
    # Step 8 — Filter network
    # -------------------------
    if should_run("filter"):
        print("[INFO] Running Step 8 — Filter network")

        grouped = pd.read_csv(
            require_file(outdir / "06_global_network_edges.tsv"),
            sep="\t"
        )
        node_base = pd.read_csv(
            require_file(outdir / "07_global_network_nodes.tsv"),
            sep="\t"
        )

        filt = grouped[grouped["n_cds_supporting_edge"] >= args.min_edge_support].copy()

        if args.min_edge_sum_cpm > 0:
            filt = filt[filt["sum_global_mean_cpm"] >= args.min_edge_sum_cpm]

        if args.top_edge_percentile > 0 and not filt.empty:
            cutoff = np.percentile(filt["sum_global_mean_cpm"], args.top_edge_percentile)
            filt = filt[filt["sum_global_mean_cpm"] >= cutoff]

        filt_nodes = sorted(set(filt["source_node"]).union(set(filt["target_node"])))
        node_filt = node_base[node_base["node_id"].isin(filt_nodes)].copy()

        filt.to_csv(outdir / "08_global_network_edges_filtered.tsv", sep="\t", index=False)
        node_filt.to_csv(outdir / "08_global_network_nodes_filtered.tsv", sep="\t", index=False)

    # -------------------------
    # Step 9 — Community detection
    # -------------------------
    if should_run("communities"):
        print("[INFO] Running Step 9 — Community detection")

        filt = pd.read_csv(
            require_file(outdir / "08_global_network_edges_filtered.tsv"),
            sep="\t"
        )
        node_filt = pd.read_csv(
            require_file(outdir / "08_global_network_nodes_filtered.tsv"),
            sep="\t"
        )

        G = nx.Graph()

        for _, r in node_filt.iterrows():
            G.add_node(
                r["node_id"],
                annotation_type=r["annotation_type"],
                annotation_id=r["annotation_id"],
                annotation_label=r["annotation_label"],
                n_cds_with_node=int(r["n_cds_with_node"]),
                sum_global_mean_cpm=float(r["sum_global_mean_cpm"]),
                median_global_mean_cpm=float(r["median_global_mean_cpm"]),
                mean_global_prevalence=float(r["mean_global_prevalence"]),
                degree=int(r["degree"]),
            )

        for _, r in filt.iterrows():
            G.add_edge(
                r["source_node"],
                r["target_node"],
                weight=float(r["sum_global_mean_cpm"]),
                edge_layer=r["edge_layer"],
                n_cds_supporting_edge=int(r["n_cds_supporting_edge"]),
                sum_global_mean_cpm=float(r["sum_global_mean_cpm"]),
                median_global_mean_cpm=float(r["median_global_mean_cpm"]),
                mean_global_prevalence=float(r["mean_global_prevalence"]),
            )

        if G.number_of_edges() > 0:
            try:
                import community as community_louvain
                partition = community_louvain.best_partition(G, weight="weight")
            except Exception:
                communities = nx.algorithms.community.greedy_modularity_communities(
                    G,
                    weight="weight"
                )
                partition = {}
                for i, comm in enumerate(communities):
                    for node in comm:
                        partition[node] = i
        else:
            partition = {}

        nx.set_node_attributes(G, partition, "community_id")

        modules = []
        for node, cid in partition.items():
            attrs = G.nodes[node]
            modules.append({
                "community_id": cid,
                "node_id": node,
                "annotation_type": attrs.get("annotation_type", ""),
                "annotation_label": attrs.get("annotation_label", ""),
            })

        modules_df = pd.DataFrame(modules)
        modules_df.to_csv(outdir / "09_global_network_modules.tsv", sep="\t", index=False)

        summaries = []
        if not modules_df.empty:
            for cid, sub in modules_df.groupby("community_id"):
                nodes_in = list(sub["node_id"])
                sub_node = node_filt[node_filt["node_id"].isin(nodes_in)]

                summaries.append({
                    "community_id": cid,
                    "n_nodes": len(nodes_in),
                    "n_pfams": int((sub["annotation_type"] == "Pfam").sum()),
                    "n_kos": int((sub["annotation_type"] == "KO").sum()),
                    "n_gos": int((sub["annotation_type"] == "GO").sum()),
                    "n_cogs": int((sub["annotation_type"] == "COG").sum()),
                    "n_ecs": int((sub["annotation_type"] == "EC").sum()),
                    "top_pfams": ";".join(
                        sub[sub["annotation_type"] == "Pfam"]["annotation_label"].head(10)
                    ),
                    "top_kos": ";".join(
                        sub[sub["annotation_type"] == "KO"]["annotation_label"].head(10)
                    ),
                    "top_go_terms": ";".join(
                        sub[sub["annotation_type"] == "GO"]["annotation_label"].head(10)
                    ),
                    "top_cog_categories": ";".join(
                        sub[sub["annotation_type"] == "COG"]["annotation_label"].head(10)
                    ),
                    "top_ecs": ";".join(
                        sub[sub["annotation_type"] == "EC"]["annotation_label"].head(10)
                    ),
                    "module_sum_global_mean_cpm": float(sub_node["sum_global_mean_cpm"].sum()),
                    "module_mean_prevalence": float(sub_node["mean_global_prevalence"].mean()),
                })

        pd.DataFrame(summaries).to_csv(
            outdir / "10_global_network_module_summary.tsv",
            sep="\t",
            index=False,
        )

    # -------------------------
    # Step 10 — Export network formats
    # -------------------------
    if should_run("export"):
        print("[INFO] Running Step 10 — Export network formats")

        filt = pd.read_csv(
            require_file(outdir / "08_global_network_edges_filtered.tsv"),
            sep="\t"
        )
        node_filt = pd.read_csv(
            require_file(outdir / "08_global_network_nodes_filtered.tsv"),
            sep="\t"
        )

        modules_path = outdir / "09_global_network_modules.tsv"
        if modules_path.is_file() and modules_path.stat().st_size > 0:
            modules_df = pd.read_csv(modules_path, sep="\t")
        else:
            modules_df = pd.DataFrame(columns=["node_id", "community_id"])

        community_map = {}
        if not modules_df.empty:
            community_map = dict(zip(modules_df["node_id"], modules_df["community_id"]))

        G = nx.Graph()

        for _, r in node_filt.iterrows():
            node_id = r["node_id"]
            G.add_node(
                node_id,
                annotation_type=r["annotation_type"],
                annotation_id=r["annotation_id"],
                annotation_label=r["annotation_label"],
                n_cds_with_node=int(r["n_cds_with_node"]),
                sum_global_mean_cpm=float(r["sum_global_mean_cpm"]),
                median_global_mean_cpm=float(r["median_global_mean_cpm"]),
                mean_global_prevalence=float(r["mean_global_prevalence"]),
                degree=int(r["degree"]),
                community_id=str(community_map.get(node_id, "NA")),
            )

        for _, r in filt.iterrows():
            G.add_edge(
                r["source_node"],
                r["target_node"],
                weight=float(r["sum_global_mean_cpm"]),
                edge_layer=r["edge_layer"],
                n_cds_supporting_edge=int(r["n_cds_supporting_edge"]),
                sum_global_mean_cpm=float(r["sum_global_mean_cpm"]),
                median_global_mean_cpm=float(r["median_global_mean_cpm"]),
                mean_global_prevalence=float(r["mean_global_prevalence"]),
            )

        nx.write_graphml(G, outdir / "global_network.graphml")
        nx.write_gexf(G, outdir / "global_network.gexf")

        print("[INFO] Export complete")
        print(f"[INFO] Nodes: {G.number_of_nodes()}")
        print(f"[INFO] Edges: {G.number_of_edges()}")

    print("[INFO] Cross-layer network workflow finished")
    print(f"[INFO] Step requested: {args.step}")
    print(f"[INFO] Output directory: {outdir}")


if __name__ == "__main__":
    main()
