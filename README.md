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
   submit_metaflye_array.sh     - submits Flye Slurm array (per-sample or global co-assembly)
   metaflye_array_task.sh       - per-array-task Flye runner (uses scratch for temp reads)
   run_smorfs_pipeline.sh       - smORFs pipeline on assemblies
   downstream_analysis.sh       - downstream analysis (Flye metrics, AMP prediction, read-back mapping)
   summarize_flye_logs.py       - parse Flye logs into a metrics table
   plot_metaflye_metrics.R      - generate per-metric boxplots from metrics TSV
   predict_amps.py              - AMP prediction using Macrel
 /envs/                         - Conda environments (created on demand)
 /logs/                         - timestamped Slurm logs
 /metadata/
   metagenome_files.txt         - SampleID ↔ FASTQ mapping (per-sample mode)
   metaflye_sampleids_*.list    - SampleID → array task mapping
   metaflye_sample_fastqs_*.tsv - SampleID ↔ FASTQ map (absolute paths, normalized)
   sample_ids.txt               - list of SampleIDs for smORFs runs
 /data/                         - input FASTQ(.gz) files (read-only; may be symlinks)
 /results/
   qc_pre_filt/                 - QC outputs on raw reads (batch-aware)
   qc_post_filt/                - QC outputs after filtering (optional)
   assembly_metaflye/           - MetaFlye assemblies (per-sample or global)
     finalize_metrics.tsv
     finalize_boxplots/
   smorfs/
   catalog_global/
 README.md
 LICENSE
 CITATION.cff
</pre>

<pre>
DESIGN PRINCIPLES
-----------------
 - Explicit separation of QC, assembly, smORFs, and downstream analysis stages
 - QC-only by default (no FASTQ files are modified)
 - Assembly is opt-in and parallelized via Slurm arrays
 - Supports both per-sample and global co-assembly modes
 - smORFs and downstream analysis are opt-in and decoupled
 - AMP prediction operates on a global non-redundant peptide catalog
 - Reproducible Conda environments (created if missing)
 - HPC-native design (Slurm + srun / sbatch)
 - Transparent logging, resumability, and provenance
 - Batch-aware QC for multi-run and longitudinal projects
 - Large temporary files are written to scratch and cleaned up automatically
</pre>

<pre>
SCRATCH USAGE
-------------
For disk safety and scalability, large temporary co-assembly FASTQ files
are written to scratch storage during MetaFlye runs and deleted automatically
after assembly completes.

Specifically:
 - Co-assembled reads are written to:
     /scratch/t.sousa/data_used/metaflye_tmp/<SampleID>/
 - These files are removed automatically after Flye finishes
 - Final assemblies, logs, and metrics are always written to the project directory
</pre>

<pre>
GLOBAL CO-ASSEMBLY MODE
----------------------
In addition to per-sample (SampleID) co-assembly, the pipeline supports
a single global metagenome assembly using all FASTQs currently present
in the data/ directory.

In global mode:
 - All libraries are merged into one co-assembly
 - A single Flye job is submitted
 - Temporary merged reads are written to scratch and deleted automatically
 - Output is written to:
     results/assembly_metaflye/<GLOBAL_ID>/
</pre>

<pre>
STEP CONTROL
------------
QC / Assembly:
  --qc-only             Run only QC (default behavior)
  --metaflye-only       Run only metagenome assembly (skip QC)
  --qc-and-metaflye     Run QC, then submit MetaFlye

MetaFlye modes:
  --global              Enable single global co-assembly across all FASTQs
  --global-id STR       SampleID name for global assembly (default: GLOBAL)

smORFs:
  --smorfs-create-env
  --smorfs-only
  --smorfs-sample STR

Downstream:
  --downstream-only
  --run-downstream
  --run-downstream-amps
  --run-downstream-map
  --downstream-full
</pre>

<pre>
ASSEMBLY OUTPUTS (MetaFlye)
--------------------------
results/assembly_metaflye/
  <SampleID>/
    assembly.fasta
    assembly_info.txt
    flye.log

Temporary (auto-deleted):
  /scratch/t.sousa/data_used/metaflye_tmp/<SampleID>/
    reads.coassembly.fastq.gz
    inputs.fastq.list
</pre>

<pre>
EXAMPLES
--------

# QC only
bash workflow/runall.sh

# QC + per-sample MetaFlye
bash workflow/runall.sh --qc-and-metaflye

# Global co-assembly (default ID = GLOBAL)
bash workflow/runall.sh --metaflye-only --global

# Global co-assembly with custom ID
bash workflow/runall.sh --metaflye-only --global --global-id UNIFORME_GLOBAL

# smORFs on one assembly
bash workflow/runall.sh --smorfs-only --smorfs-sample TS-0500
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
