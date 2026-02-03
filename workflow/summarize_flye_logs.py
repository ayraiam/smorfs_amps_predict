#!/usr/bin/env python3
import re
import sys
from pathlib import Path
import pandas as pd

def parse_flye_log(log_path: Path) -> dict:
    """
    Best-effort parser for Flye/MetaFlye flye.log.
    Scans lines and matches known patterns; keeps last match per metric.
    """
    text = log_path.read_text(errors="replace").splitlines()

    out = {
        "site": log_path.parent.name,
        "total_reads_bp": None,
        "total_assembly_mb": None,
        "contigs": None,
        "n50_bp": None,
        "largest_bp": None,
        "mean_cov": None,
        "flye_log": str(log_path),
    }

    patterns = {
        "total_reads_bp": [
            re.compile(r"Total read length:\s*([0-9]+)"),
            re.compile(r"Total reads length:\s*([0-9]+)"),
        ],
        "total_assembly_bp": [
            re.compile(r"Total length:\s*([0-9]+)"),
            re.compile(r"Total length of contigs:\s*([0-9]+)"),
            re.compile(r"Total assembled length:\s*([0-9]+)"),
        ],
        "contigs": [
            # Flye often prints something like:
            #   Fragments: <N>
            # or
            #   Fragments <N>
            re.compile(r"^\s*Fragments\s*[:=]\s*([0-9]+)\s*$"),
            re.compile(r"^\s*Fragments\s+([0-9]+)\s*$"),
        ],
        "n50_bp": [
            # Example:
            #   Fragments N50: <N>
            re.compile(r"^\s*Fragments\s+N50\s*[:=]\s*([0-9]+)\s*$"),
        ],
        "largest_bp": [
            # Example:
            #   Largest frg: <N>
            re.compile(r"^\s*Largest\s+frg\s*[:=]\s*([0-9]+)\s*$"),
        ],

        "mean_cov": [
            re.compile(r"Mean coverage:\s*([0-9]+(?:\.[0-9]+)?)"),
            re.compile(r"Mean cov(?:erage)?:\s*([0-9]+(?:\.[0-9]+)?)"),
        ],
    }

    assembly_bp = None

    for line in text:
        line = line.strip()

        for rx in patterns["total_reads_bp"]:
            m = rx.search(line)
            if m:
                out["total_reads_bp"] = int(m.group(1))

        for rx in patterns["total_assembly_bp"]:
            m = rx.search(line)
            if m:
                assembly_bp = int(m.group(1))

        for rx in patterns["contigs"]:
            m = rx.search(line)
            if m:
                out["contigs"] = int(float(m.group(1)))

        for rx in patterns["n50_bp"]:
            m = rx.search(line)
            if m:
                out["n50_bp"] = int(m.group(1))

        for rx in patterns["largest_bp"]:
            m = rx.search(line)
            if m:
                out["largest_bp"] = int(m.group(1))

        for rx in patterns["mean_cov"]:
            m = rx.search(line)
            if m:
                out["mean_cov"] = float(m.group(1))

    if assembly_bp is not None:
        out["total_assembly_mb"] = assembly_bp / 1e6

    return out

def usage():
    print(
        "Usage: summarize_flye_logs.py <BASE_RESULTS_DIR> [--out <tsv>] [--no-plots]\n"
        "\n"
        "Parses Flye/MetaFlye flye.log files under:\n"
        "  <BASE_RESULTS_DIR>/assembly_metaflye/*/flye.log\n"
        "\n"
        "Writes a TSV table (one row per SampleID/site).\n"
    )

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        usage()
        sys.exit(0 if len(sys.argv) >= 2 else 1)

    base_results = Path(sys.argv[1]).resolve()

    # Default output location
    out_tsv = base_results / "assembly_metaflye" / "finalize_metrics.tsv"

    # Parse args (keep --no-plots for backward compatibility)
    argv = sys.argv[2:]
    if "--out" in argv:
        i = argv.index("--out")
        try:
            out_tsv = Path(argv[i + 1]).resolve()
        except IndexError:
            raise SystemExit("ERROR: --out requires a path")

    logs = sorted(base_results.glob("assembly_metaflye/*/flye.log"))
    if not logs:
        raise SystemExit(f"No flye.log found under: {base_results}/assembly_metaflye/*/flye.log")

    rows = [parse_flye_log(p) for p in logs]
    df = pd.DataFrame(rows)

    # Stable column order
    cols = ["site", "total_reads_bp", "total_assembly_mb", "contigs", "n50_bp", "largest_bp", "mean_cov", "flye_log"]
    df = df.reindex(columns=[c for c in cols if c in df.columns])

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"Wrote: {out_tsv}")
    print(f"Rows:  {df.shape[0]}")

if __name__ == "__main__":
    main()
