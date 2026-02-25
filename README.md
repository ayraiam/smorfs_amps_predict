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
smorfs_amps_predict is a lightweight, reproducible, HPC-native pipeline for
Oxford Nanopore shotgun metagenomic data, providing a modular framework for:

  1) Read-level quality control and diagnostics (safe, read-only by default)
  2) Scalable metagenome co-assembly using Flye (--meta)
     • Per-sample (SampleID) co-assemblies
     • Single global co-assembly across all libraries
  3) On-demand smORF prediction (bacterial + fungal)
  4) AMP prediction via Macrel (two operational modes)
  5) Bacterial smORF refinement based on genomic context
  6) Downstream assembly metrics + visualization

All major stages are explicitly decoupled:
  - QC once → reuse
  - Assembly without repeating QC
  - smORFs independently
  - AMP prediction independently
  - Refinement independently
  - Downstream analysis independently

This design enables reproducibility, re-entrancy, and HPC-scalable workflows.
</pre>

---

<pre>
STRUCTURE
---------
 /workflow/
   runall.sh                               - Main Slurm launcher
   run_libsQC.sh                           - QC logic (FastQC, NanoPlot, NanoStat, SeqKit)
   submit_metaflye_array.sh                - Flye Slurm array submission
   metaflye_array_task.sh                  - Per-array Flye runner (scratch-aware)
   run_smorfs_pipeline.sh                  - smORF discovery pipeline
   smorfs_job.sh                           - Slurm wrapper for smORFs
   run_refine_annot_smorf_bacs.sh          - Bacterial refinement stage
   refine_bacs_job.sh                      - Slurm wrapper for refinement
   refine_bacs.py                          - Genomic-context refinement logic
   downstream_analysis.sh                  - Flye metrics + AMP utilities
   summarize_flye_logs.py                  - Parse Flye logs into metrics table
   plot_metaflye_metrics.R                 - Boxplots from metrics TSV
   predict_amps.py                         - Macrel-based AMP prediction
   attach_macrel_to_predicted_smorfs.py    - Safe Macrel merge (ID-validated)

 /envs/                                    - Conda environments (auto-created)
 /logs/                                    - Timestamped Slurm logs
 /metadata/
   metagenome_files.txt
   metaflye_sampleids_*.list
   metaflye_sample_fastqs_*.tsv
   sample_ids.txt
 /data/                                    - Input FASTQ(.gz) files (read-only)
 /results/
   qc_pre_filt/
   qc_post_filt/
   assembly_metaflye/
   smorfs/
   catalog_global/
 README.md
 LICENSE
 CITATION.cff
</pre>

---

<pre>
DESIGN PRINCIPLES
-----------------
 - Explicit stage separation (QC → Assembly → smORFs → AMP → Refinement → Downstream)
 - QC-only by default (no FASTQs modified)
 - Assembly opt-in, parallelized via Slurm arrays
 - Supports per-sample and global co-assembly
 - smORFs and refinement fully decoupled from assembly
 - Deterministic ID-safe Macrel merging
 - Hard failure on sequence mismatches (no silent corruption)
 - Automatic Conda environment creation
 - Slurm-native orchestration (srun + sbatch)
 - Transparent logging, resumability, provenance tracking
 - Scratch-aware temporary file handling
 - Batch-aware QC for multi-run projects
</pre>

---

<pre>
BACTERIAL smORF REFINEMENT
--------------------------
After smORF prediction (and optional Macrel AMP annotation),
the pipeline supports a bacterial refinement stage that:

  • Integrates Prodigal GFF annotations
  • Integrates bacterial contig FASTA
  • Builds genomic-context dictionaries
  • Evaluates structural context of predicted smORFs
  • Produces refined TSV outputs

Current refinement framework includes fields for:

  flag_edge
  dist_left / dist_right
  flag_embedded
  flag_overlap_fraction
  host_cds_id
  cluster_id / cluster_size
  cluster_env_count
  flag_cluster_recurrent
  flag_cross_environment
  flag_too_short
  confidence_tier

Refinement is designed to classify predictions into confidence tiers
based on genomic context and recurrence, enabling structured downstream filtering.

Refinement is independent and can be run:

  • After smORFs
  • After Macrel annotation
  • On existing catalogs
  • For one sample or multiple samples

If refinement is launched in the same run as smORFs,
Slurm dependency is automatically enforced (afterok).
</pre>

---

<pre>
MACREL AMP MODES
----------------
Two AMP prediction strategies are supported:

1) Global NR Catalog Mode
   - Concatenate peptide catalogs
   - Cluster with mmseqs
   - Run Macrel on non-redundant representatives
   - Output:
       results/catalog_global/global_nr_peptides.faa
       results/catalog_global/amp_predictions.tsv

2) Direct Per-Sample Annotation Mode
   - Input:
       results/smorfs/<SampleID>/catalog/predicted_smorfs.tsv
   - Builds FASTA using stable ID column
   - Runs Macrel
   - Validates sequence identity before merging
   - Output:
       predicted_smorfs.with_macrel.tsv

Macrel fields added:
  amp_pred
  amp_prob
  hemo_pred
  hemo_prob
  macrel_class
</pre>

---

<pre>
STEP CONTROL
------------

QC / Assembly:
  --qc-only
  --metaflye-only
  --qc-and-metaflye

MetaFlye:
  --global
  --global-id STR

smORFs:
  --smorfs-only
  --smorfs-sample STR
  --smorfs-samples FILE
  --smorfs-create-env

Refinement (Bacterial):
  --refine-bacs-only
  --run-refine-bacs
  --refine-bacs-sample STR
  --refine-bacs-samples FILE
  --refine-bacs-create-env

Downstream:
  --downstream-only
  --run-downstream
  --run-downstream-amps
  --run-downstream-map
  --downstream-full

Macrel Attach Mode:
  --macrel-attach-only
  --predicted-smorfs PATH
  --id-col STR
  --seq-col STR
  --macrel-attach-out PATH
</pre>

---

<pre>
SCRATCH USAGE
-------------
Large temporary co-assembly FASTQ files are written to scratch:

  /scratch/t.sousa/data_used/metaflye_tmp/<SampleID>/

These are deleted automatically after Flye completes.

Final assemblies and outputs are always written inside the project directory.
</pre>

---

<pre>
EXAMPLES
--------

# 1) QC only (default behavior)
bash workflow/runall.sh

# 2) QC + per-sample MetaFlye assembly
bash workflow/runall.sh --qc-and-metaflye

# 3) Global co-assembly across all FASTQs
bash workflow/runall.sh --metaflye-only --global

# 4) Run smORFs on one assembly
bash workflow/runall.sh \
  --smorfs-only \
  --smorfs-sample TS-0500

# 5) Annotate smORFs with Macrel (direct per-sample mode)
bash workflow/runall.sh \
  --macrel-attach-only \
  --predicted-smorfs results/smorfs/TS-0500/catalog/predicted_smorfs.tsv \
  --id-col feature_id \
  --seq-col aa_seq

# 6) Run bacterial refinement (Step 1: edge flagging)
bash workflow/runall.sh \
  --refine-bacs-only \
  --refine-bacs-sample TS-0500

# 7) Run smORFs + refinement in one invocation (dependency enforced)
bash workflow/runall.sh \
  --smorfs-only \
  --smorfs-sample TS-0500 \
  --run-refine-bacs

# 8) Run refinement for multiple samples
bash workflow/runall.sh \
  --refine-bacs-only \
  --refine-bacs-samples metadata/sample_ids.txt

# 9) Full structured workflow (assembly → smORFs → refinement)
bash workflow/runall.sh \
  --metaflye-only --global \
  --run-smorfs \
  --run-refine-bacs

# 10) Downstream assembly metrics + AMP analysis
bash workflow/runall.sh --downstream-full
</pre>

---

<pre>
CITATION
--------
Lobo, I. (2026).
smorfs_amps_predict: A reproducible, HPC-native pipeline for
quality control, metagenome assembly, smORF discovery,
AMP prediction, and genomic-context refinement
of Nanopore shotgun data.
AY:RΔ — data and discovery in flow.
</pre>

<p align="center"><sub>© 2026 AY:RΔ — data and discovery in flow</sub></p>
