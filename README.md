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
smorfs_amps_predict is a lightweight, reproducible pipeline for
Oxford Nanopore shotgun metagenomic data, providing a robust backbone
for downstream smORF and antimicrobial peptide (AMP) discovery.

The pipeline combines:
  1) Read-level quality control and diagnostics (safe, read-only by default)
  2) Scalable metagenome assembly using Flye (--meta), parallelized via Slurm arrays

QC and assembly are explicitly decoupled, allowing users to re-run
assembly without repeating QC, resume partial runs, and efficiently
scale large metagenome collections on HPC systems.

QC diagnostics are batch-aware, enabling direct comparison of
multiple sequencing batches without overwriting results.
</pre>

<pre>
STRUCTURE
---------
 /workflow/
   runall.sh                    - main Slurm launcher (step control)
   run_libsQC.sh                - QC logic (FastQC, NanoPlot, NanoStat, SeqKit)
   submit_metaflye_array.sh     - submits Flye Slurm array (1 SampleID per task; co-assembly)
   metaflye_array_task.sh       - per-array-task Flye runner (co-assembly)
 /envs/                         - Conda environments (created on demand)
 /logs/                         - timestamped Slurm logs
 /metadata/
   metagenome_files.txt         - SampleID ↔ FASTQ mapping (user-provided)
   metaflye_sampleids_*.list    - SampleID → array task mapping
   metaflye_sample_fastqs_*.tsv - SampleID ↔ FASTQ map (absolute paths, normalized)
 /data/                         - input FASTQ(.gz) files (read-only)
 /results/
   qc_pre_filt/                 - QC outputs on raw reads (batch-aware)
   qc_post_filt/                - QC outputs after filtering (optional; batch-aware)
   assembly_metaflye/           - per-sample Flye metagenome assemblies
 README.md                      - this file
 LICENSE                        - project license (MIT)
 CITATION.cff                   - citation metadata
</pre>

<pre>
DEPENDENCIES
------------
 - Linux
 - Bash ≥ 4
 - Conda or Mamba
 - Slurm (or compatible scheduler)

Automatically installed in pipeline environments:
 - FastQC
 - MultiQC
 - NanoPlot
 - NanoStat
 - SeqKit
 - Cutadapt
 - NanoFilt
 - Porechop
 - Flye
</pre>

<pre>
DESIGN PRINCIPLES
-----------------
 - Explicit separation of QC and assembly stages
 - QC-only by default (no FASTQ files are modified)
 - Assembly is opt-in and fully parallelized
 - Co-assembly is performed per biological sample (SampleID)
 - Reproducible Conda environments (created if missing)
 - HPC-native design (Slurm + srun / sbatch arrays)
 - Transparent logging, resumability, and provenance
 - Batch-aware QC for multi-run and longitudinal projects

</pre>

<pre>
BATCH AWARE QC
--------------
The QC stage supports explicit batch labeling via a user-defined
batch identifier.

This enables:
 - Direct comparison of multiple sequencing batches
 - Preservation of independent NanoPlot / NanoStat outputs
 - Safe re-analysis without overwriting previous QC results

Batch labeling affects:
 - NanoPlot outputs
 - NanoStat summaries
 - SeqKit read statistics

No FASTQ files are ever modified.
</pre>

<pre>
USAGE
-----
Place FASTQ or FASTQ.GZ files in the `data/` directory.

Provide a mapping file at:
  metadata/metagenome_files.txt

This file must associate each FASTQ file with a SampleID.
It can be tab-, comma-, or space-delimited, with or without a header.

Example:
  SampleID    FASTQ_Filename
  S01         lib1.fastq.gz
  S01         lib2.fastq.gz
  S02         lib3.fastq.gz

Run the main launcher:

  bash workflow/runall.sh [options]

Optionally, specify a batch identifier to keep QC outputs from
different sequencing runs separate.

By default:
 - Only the QC stage is executed
 - FASTQ files are never modified
 - No assembly is performed unless explicitly requested
</pre>

<pre>
STEP CONTROL
------------
The pipeline supports explicit step selection:

  --qc-only             Run only QC (default behavior)
  --metaflye-only       Run only metagenome assembly (skip QC)
  --qc-and-metaflye     Run QC, then submit Flye assemblies

This allows safe re-use of QC results and re-running assembly
without repeating earlier steps.
</pre>

<pre>
MAIN OPTIONS
------------
Resources:
  --partition STR       Slurm partition (default: short)
  --time HH:MM:SS       Walltime per job/task
  --cpus INT            CPUs per task
  --mem STR             Memory per task
  --wd PATH             Working directory

QC diagnostics:
  --run-porechop        Run porechop adapter/barcode diagnostics (OFF by default)

Filtering (used only if --run-filtering is enabled):
  --run-filtering       Enable trimming/filtering stage
  --min-q INT           Mean read Q cutoff (default: 10)
  --min-len INT         Minimum read length (default: 500 bp)
  --max-len INT         Maximum read length (0 = disabled)

  --no-adapter-trim     Skip adapter trimming
  --no-barcode-trim     Skip barcode trimming
  --demux               Enable demultiplexing
  --no-poly-trim        Skip poly-A/T trimming
  --no-filter           Skip NanoFilt Q/length filtering

Batch control:
  --batch-id STR        Batch identifier for QC outputs
                        (default: batch1)

</pre>

<pre>
QC OUTPUTS
----------
results/qc_pre_filt/
  fastqc/                         - per-file FastQC reports
  multiqc/                        - aggregated MultiQC report

  nanoplot/<batch>_raw/           - read length & quality distributions
  nanostat/<batch>_raw/           - per-file NanoStat summaries
  summary/
    seqkit_stats_<batch>_raw.tsv  - SeqKit read statistics

  adapter_barcode_checks/         - porechop diagnostics (if enabled)

results/qc_post_filt/
  nanoplot/<batch>_clean/
  nanostat/<batch>_clean/
  summary/
    seqkit_stats_<batch>_clean.tsv

Where <batch> corresponds to the value provided via --batch-id.
</pre>

<pre>
ASSEMBLY OUTPUTS (MetaFlye)
--------------------------
results/assembly_metaflye/
  <SampleID>/
    reads.coassembly.fastq.gz  - concatenated reads used for co-assembly
    inputs.fastq.list          - FASTQs included in the co-assembly
    assembly.fasta             - assembled contigs
    assembly_info.txt          - contig stats and coverage
    flye.log                   - Flye internal log
  ...

Reads are co-assembled per SampleID: all FASTQs mapped to the same
SampleID (via metadata/metagenome_files.txt) are concatenated and
assembled using Flye --meta.

Assemblies are skipped automatically if output already exists,
allowing safe re-runs and recovery from partial failures.
</pre>

<pre>
PARALLELISM MODEL
-----------------
 - One Slurm array task per SampleID
 - Each task concatenates that SampleID’s FASTQs into a single
   co-assembly input and runs Flye
 - Each task receives its own CPUs, memory, and walltime
 - Conda environment creation is lock-protected and occurs only once
 - Designed for large-scale metagenome collections
</pre>

<pre>
LOGGING & REPRODUCIBILITY
------------------------
 - logs/qc_*.out / *.err
     QC Slurm stdout/stderr

 - logs/metaflye_array_*_<job>_<task>.out/.err
     Per-sample assembly logs

 - logs/command_*.txt
     Full pipeline invocation

 - metadata/metagenome_files.txt
     User-provided SampleID ↔ FASTQ mapping

 - metadata/metaflye_sampleids_*.list
     SampleID → Slurm array task mapping

 - metadata/metaflye_sample_fastqs_*.tsv
     Normalized SampleID ↔ FASTQ map (absolute paths)

 - envs/
     Conda environments (created on demand)
</pre>

<pre>
EXAMPLES
--------

# 1) Default QC-only run (safe, read-only)
bash workflow/runall.sh

# 2) QC-only run with explicit batch labeling
bash workflow/runall.sh --batch-id batch1

# 3) QC + porechop diagnostics (batch-aware)
bash workflow/runall.sh --run-porechop --batch-id batch2

# 4) MetaFlye co-assembly only (skip QC)
bash workflow/runall.sh --metaflye-only --cpus 16 --mem 64G --time 12:00:00

# 5) QC followed by MetaFlye co-assembly
bash workflow/runall.sh --qc-and-metaflye --batch-id batch1 --cpus 16 --mem 64G

# 6) Advanced: QC + filtering + co-assembly
bash workflow/runall.sh --qc-and-metaflye --run-filtering \
  --batch-id batch2 --min-q 12 --min-len 1000

</pre>

<pre>
NOTES ON PORECHOP
----------------
Porechop is used ONLY as a diagnostic tool when --run-porechop is enabled.
It does NOT modify FASTQ files.

This is intentional, as modern ONT basecalling software (e.g. MinKNOW)
already removes adapters in most workflows. The porechop step helps
identify residual or unexpected adapter signals.
</pre>

<pre>
CITATION
--------
If you use this pipeline, please cite:

Lobo, I. (2026).
smorfs_amps_predict: A reproducible pipeline for quality control
and metagenome assembly of Nanopore shotgun data.
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

