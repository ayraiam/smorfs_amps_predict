<!-- LOGO -->
<p align="center">
  <img src="https://github.com/ayrdelta/.github/blob/main/profile/logo.png" alt="AY:RΔ Logo" width="120"/>
</p>

<h1 align="center" style="font-weight: normal;">smorfs_amps_predict</h1>

<p align="center">
  <code>data and discovery in flow</code><br/>
  <a href="mailto:ayrabioinf@gmail.com">ayrabioinf@gmail.com</a> · 
  <a href="https://www.linkedin.com/company/aryaiam">LinkedIn</a>
</p>

---

<pre>
ABOUT
-----
smorfs_amps_predict is a lightweight, reproducible quality-control pipeline
for Oxford Nanopore shotgun metagenomic FASTQ data.

The pipeline focuses on read-level QC and diagnostics without modifying
the original FASTQ files by default, making it safe for exploratory
analysis, troubleshooting, and reporting.

It was designed for large-scale environmental and clinical metagenome
projects running on HPC systems with Slurm.
</pre>

<pre>
STRUCTURE
---------
 /workflow/
   runall.sh          - main Slurm launcher
   run_libsQC.sh      - QC logic (FastQC, NanoPlot, NanoStat, SeqKit)
 /envs/               - exported Conda environments
 /logs/               - timestamped logs and timing information
 /metadata/           - reserved for sample / run metadata
 /data/               - input FASTQ(.gz) files (read-only)
 /results/
   qc_pre_filt/       - QC outputs on raw reads
   qc_post_filt/      - QC outputs after filtering (optional)
 README.md            - this file
 LICENSE              - project license (MIT)
 CITATION.cff         - citation metadata
</pre>

<pre>
DEPENDENCIES
------------
 - Linux / macOS
 - Bash ≥ 4
 - Conda or Mamba
 - Slurm (or compatible scheduler)

Automatically installed in the pipeline environment:
 - FastQC
 - MultiQC
 - NanoPlot
 - NanoStat
 - SeqKit
 - Cutadapt
 - NanoFilt
 - Porechop
</pre>

<pre>
DESIGN PRINCIPLES
-----------------
 - QC-only by default (no FASTQ files are modified)
 - Explicit opt-in for any potentially destructive operation
 - Reproducible environments (exported YAML)
 - HPC-friendly (Slurm + srun)
 - Transparent logging and timing
</pre>

<pre>
USAGE
-----
Place your FASTQ or FASTQ.GZ files in the `data/` directory and run:

  bash workflow/runall.sh [options]

By default, the pipeline:
 - Creates a Conda environment (metaQC)
 - Runs read-only QC
 - Produces summary plots and tables
 - Does NOT trim, filter, or modify reads
</pre>

<pre>
MAIN OPTIONS
------------
Resources:
  --partition STR       Slurm partition (default: short)
  --time HH:MM:SS       Walltime (default: 04:00:00)
  --cpus INT            CPUs per task (default: 8)
  --mem STR             Memory (default: 32G)
  --wd PATH             Working directory (default: current directory)

Mode:
  --run-filtering       Enable filtering/cleanup stage (OFF by default)

QC checks:
  --run-porechop        Run porechop adapter/barcode diagnostic checks
                        (OFF by default)

Filtering parameters (used only if --run-filtering is enabled):
  --min-q INT           Mean read Q cutoff (default: 10)
  --min-len INT         Minimum read length in bp (default: 500)
  --max-len INT         Maximum read length (0 = disabled)

  --no-adapter-trim     Skip adapter trimming
  --no-barcode-trim     Skip barcode trimming
  --demux               Enable demultiplexing
  --no-poly-trim        Skip poly-A/T trimming
  --no-filter           Skip NanoFilt Q/length filtering
</pre>

<pre>
QC OUTPUTS
----------
results/qc_pre_filt/
  fastqc/                     - per-file FastQC reports
  multiqc/                    - aggregated MultiQC report
  nanoplot/raw/               - read length & quality distributions
  nanostat/raw/               - per-sample NanoStat summaries
  summary/                    - SeqKit read statistics
  adapter_barcode_checks/     - porechop diagnostics (if enabled)

results/qc_post_filt/
  (same structure, only if filtering is enabled)
</pre>

<pre>
LOGGING & REPRODUCIBILITY
------------------------
 - logs/qc_*.out / *.err      Slurm stdout/stderr
 - logs/command_*.txt         Full invocation record
 - logs/.timing.tsv           Per-step runtime profiling
 - envs/metaQC.yml            Exported Conda environment
</pre>

<pre>
EXAMPLES
--------

# 1) Default QC-only run (safe, read-only)
bash workflow/runall.sh

# 2) QC + porechop adapter/barcode diagnostics
bash workflow/runall.sh --run-porechop

# 3) Increase resources
bash workflow/runall.sh --time 06:00:00 --cpus 16 --mem 64G

# 4) QC-only on a specific working directory
bash workflow/runall.sh --wd /path/to/project

# 5) Enable filtering + post-QC (advanced / optional)
bash workflow/runall.sh --run-filtering --min-q 12 --min-len 1000
</pre>

<pre>
NOTES ON PORECHOP
----------------
Porechop is used here ONLY as a diagnostic tool when --run-porechop
is enabled. It does NOT trim reads or modify FASTQs.

This is intentional, as modern ONT basecalling software (e.g. MinKNOW)
already removes adapters in most workflows. The porechop step helps
identify residual or unexpected adapter signals.
</pre>

<pre>
CITATION
--------
If you use this pipeline, please cite:

Lobo, I. (2026).
smorfs_amps_predict: A reproducible quality-control pipeline
for Nanopore shotgun metagenomics.
AY:RΔ — data and discovery in flow.

See CITATION.cff for full citation metadata.
</pre>

<pre>
CONTACT
-------
ayrabioinf@gmail.com
https://www.linkedin.com/company/aryaiam
</pre>

<p align="center"><sub>© 2026 AY:RΔ — data and discovery in flow</sub></p>
