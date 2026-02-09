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
  2) Scalable metagenome co-assembly using Flye (--meta), supporting both:
     - per-sample (SampleID) co-assemblies
     - single global co-assembly across all libraries
     Parallelized via Slurm arrays
  3) Optional, on-demand smORF prediction (bacterial + fungal) on existing assemblies
  4) Optional downstream analysis of MetaFlye assemblies (flye.log → metrics TSV + boxplots)

QC, assembly, smORF prediction, and downstream analysis are explicitly decoupled:
  - Run QC once, reuse outputs
  - Re-run assembly without repeating QC
  - Run smORFs later (or re-run) without repeating QC/assembly
  - Run downstream analysis after assemblies are available (no FASTQs required)

QC diagnostics are batch-aware, enabling direct comparison of
multiple sequencing batches without overwriting results.
</pre>

<pre>
STRUCTURE
---------
 /workflow/
   runall.sh                    - main Slurm launcher (QC + assembly; optional smORFs submission)
   run_libsQC.sh                - QC logic (FastQC, NanoPlot, NanoStat, SeqKit)
   submit_metaflye_array.sh     - submits Flye Slurm array (1 SampleID per task; co-assembly)
   metaflye_array_task.sh       - per-array-task Flye runner (co-assembly; uses scratch for temp reads)
   run_smorfs_pipeline.sh       - smORFs pipeline on assemblies (Tiara → Prodigal/SmORFinder or funannotate → MMseqs2)
   downstream_analysis.sh       - downstream analysis (Flye metrics, AMP prediction, read-back mapping)
   summarize_flye_logs.py       - parse Flye logs into a metrics table
   plot_metaflye_metrics.R      - generate per-metric boxplots from metrics TSV
   predict_amps.py              - AMP prediction using Macrel (global NR peptides)
 /envs/                         - Conda environments (created on demand)
 /logs/                         - timestamped Slurm logs
 /metadata/
   metagenome_files.txt         - SampleID ↔ FASTQ mapping (user-provided)
   metaflye_sampleids_*.list    - SampleID → array task mapping
   metaflye_sample_fastqs_*.tsv - SampleID ↔ FASTQ map (absolute paths, normalized)
   sample_ids.txt               - list of SampleIDs for smORFs runs (one per line)
 /data/                         - input FASTQ(.gz) files (read-only)
 /results/
   qc_pre_filt/                 - QC outputs on raw reads (batch-aware)
   qc_post_filt/                - QC outputs after filtering (optional; batch-aware)
   assembly_metaflye/           - per-sample Flye metagenome assemblies
     finalize_metrics.tsv       - downstream: aggregated Flye assembly metrics
     finalize_boxplots/         - downstream: boxplots summarizing Flye metrics
   smorfs/                      - per-sample smORF outputs + global summary TSVs
   catalog_global/              - global NR peptides, AMP predictions, mapping outputs
 README.md
 LICENSE
 CITATION.cff
</pre>

<pre>
DEPENDENCIES
------------
 - Linux
 - Bash ≥ 4
 - Conda or Mamba
 - Slurm (or compatible scheduler)

QC + assembly environments (created automatically as needed) install:
 - FastQC
 - MultiQC
 - NanoPlot
 - NanoStat
 - SeqKit
 - Cutadapt
 - NanoFilt
 - Porechop (diagnostic only)
 - Flye (metagenome mode; Flye --meta)

smORFs environment (created on demand) installs:
 - Prodigal (meta mode)
 - Tiara (contig classification; euk vs prok)
 - SmORFinder (bacterial smORF enrichment/annotation)
 - funannotate (eukaryotic gene prediction on fungal contigs)
 - MMseqs2 (non-redundant peptide catalog)
 - SeqKit (FASTA utilities)

Downstream analysis environment (created on first run of downstream_analysis.sh) installs:
 - Python ≥ 3.10
 - pandas
 - matplotlib
 - R ≥ 4.3
 - ggplot2
 - readr
 - dplyr
 - ggbeeswarm
 - MMseqs2
 - minimap2
 - samtools
 - bedtools
 - seqkit
 - macrel
</pre>

<pre>
DESIGN PRINCIPLES
-----------------
 - Explicit separation of QC, assembly, smORFs, and downstream analysis stages
 - QC-only by default (no FASTQ files are modified)
 - Assembly is opt-in and fully parallelized (Slurm arrays)
 - Co-assembly is performed per biological sample (SampleID)
 - smORFs are opt-in and run only when explicitly requested
 - Downstream analysis is opt-in and can run independently
 - AMP prediction is performed on a global non-redundant peptide catalog
 - Read-back mapping is replicate-aware and decoupled from assembly
 - Reproducible Conda environments (created if missing)
 - HPC-native design (Slurm + srun / sbatch arrays)
 - Transparent logging, resumability, and provenance
 - Batch-aware QC for multi-run and longitudinal projects
 - Large temporary files are written to scratch storage and cleaned up automatically
</pre>

<pre>
SCRATCH USAGE
------------------------
For disk safety and scalability, large temporary co-assembly FASTQ files
are written to scratch storage during MetaFlye runs and deleted automatically
after assembly completes.

Specifically:
 - Co-assembled reads (reads.coassembly.fastq.gz) are written to:
     /scratch/t.sousa/data_used/metaflye_tmp/<SampleID>/
 - These files are removed automatically after Flye finishes
 - Final assemblies, logs, and metrics are always written to the project directory

This design prevents project disk exhaustion while preserving full reproducibility.
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

By default:
 - Only the QC stage is executed
 - FASTQ files are never modified
 - No assembly is performed unless explicitly requested
 - No smORFs are run unless explicitly requested
</pre>

<pre>
GLOBAL CO-ASSEMBLY MODE
----------------------
In addition to per-sample (SampleID) co-assembly, the pipeline supports
a single global metagenome assembly using all FASTQs currently present
in the data/ directory.

In global mode:
 - All libraries are merged into one co-assembly
 - A single Flye run is submitted
 - Temporary merged reads are written to scratch and deleted after completion
 - Output is written to:
     results/assembly_metaflye/GLOBAL/
</pre>

<pre>
STEP CONTROL
------------
QC / Assembly:
  --qc-only             Run only QC (default behavior)
  --metaflye-only       Run only metagenome assembly (skip QC)
  --qc-and-metaflye     Run QC, then submit Flye assemblies

smORFs (on existing assemblies):
  --smorfs-create-env   Create smORFs environment and exit
  --smorfs-only         Submit smORFs job (assumes MetaFlye outputs exist)
  --smorfs-sample STR   Run smORFs for ONE SampleID only
  --run-smorfs          Compatibility flag (legacy)

Downstream:
  --downstream-only          Run downstream analysis only
  --run-downstream           Run downstream analysis in addition to selected steps
  --run-downstream-amps      Enable AMP prediction (Macrel)
  --run-downstream-map       Enable read-back mapping for abundance estimation
  --downstream-full          Run downstream analysis + AMP prediction + mapping
  --metrics-env STR          Env name for downstream analysis
</pre>

<pre>
ASSEMBLY OUTPUTS (MetaFlye)
--------------------------
results/assembly_metaflye/
  <SampleID>/
    assembly.fasta
    assembly_info.txt
    flye.log

Temporary (auto-deleted) during assembly:
  /scratch/t.sousa/data_used/metaflye_tmp/<SampleID>/
    reads.coassembly.fastq.gz
    inputs.fastq.list

</pre>

<pre>
DOWNSTREAM OUTPUTS (MetaFlye logs)
---------------------------------
Inputs:
  results/assembly_metaflye/<SampleID>/flye.log

Outputs:
  results/assembly_metaflye/finalize_metrics.tsv
  results/assembly_metaflye/finalize_boxplots/boxplot_*.png
</pre>

<pre>
DOWNSTREAM OUTPUTS (EXTENDED: AMP + ABUNDANCE)
---------------------------------------------
Additional downstream steps operate on existing assemblies and smORF outputs.

Global non-redundant peptide catalog:
  results/catalog_global/global_nr_peptides.faa

AMP prediction (Macrel):
  results/catalog_global/amp_predictions.tsv
  results/catalog_global/macrel_out_peptides/

Read-back mapping (per replicate):
  results/catalog_global/mapping/
    <replicate>.idxstats.tsv
    <replicate>.mean_depth_per_contig.tsv

Mapping is performed per replicate, enabling aggregation
at the SampleID or environment level downstream.
</pre>

<pre>
EXAMPLES
--------

# Default QC-only run
bash workflow/runall.sh

# QC + MetaFlye
bash workflow/runall.sh --qc-and-metaflye

# smORFs on one assembly
bash workflow/runall.sh --smorfs-only --smorfs-sample TS-0500

# Downstream metrics only
bash workflow/runall.sh --downstream-only

# Downstream + AMP prediction
bash workflow/runall.sh --downstream-only --run-downstream-amps

# Full downstream (metrics + AMP + mapping)
bash workflow/runall.sh --downstream-only --downstream-full
</pre>

<pre>
CITATION
--------
Lobo, I. (2026).
smorfs_amps_predict: A reproducible pipeline for quality control,
metagenome assembly, and on-demand smORF and AMP discovery of
Nanopore shotgun data.
AY:RΔ — data and discovery in flow.
</pre>

<p align="center"><sub>© 2026 AY:RΔ — data and discovery in flow</sub></p>





