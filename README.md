<!-- LOGO -->
<p align="center">
  <img src="assets/logo-horizontal.png" alt="AY:RΔ Logo" width="160"/>
</p>

<h1 align="center" style="font-weight: normal;">smorfs_amps_predict</h1>

<p align="center">
  <img src="https://img.shields.io/badge/pipeline-HPC--native-green"/>
  <img src="https://img.shields.io/badge/data-Nanopore-orange"/>
  <img src="https://img.shields.io/badge/target-smORFs%20%2B%20AMPs-purple"/>
  <img src="https://img.shields.io/badge/license-MIT-blue"/>
</p>

<p align="center">
  <code>data and discovery in flow</code><br/>
  <a href="mailto:ayrabioinf@gmail.com">ayrabioinf@gmail.com</a> · 
  <a href="https://www.linkedin.com/company/aryaiam">LinkedIn</a>
</p>
---

<pre>
ABOUT
-----
smorfs_amps_predict is a reproducible, HPC-native pipeline for
Oxford Nanopore shotgun metagenomic data, providing a modular framework for:

  1) Read-level quality control and diagnostics (safe, read-only by default)
  2) Scalable metagenome co-assembly using Flye (--meta)
     • Per-sample (SampleID) co-assemblies
     • Single global co-assembly across all libraries
  3) smORF and short-protein discovery
   • Bacterial ORFs via Prodigal
   • smORF candidates via SmORFinder
   • Eukaryotic/fungal proteins via MetaEuk
  4) AMP prediction via Macrel
   • Supports bacterial/archaeal/prokarya and fungal peptides
  5) smORF refinement based on genomic context
   • bacterial refinement (Prodigal context)
   • fungal/euk refinement (MetaEuk context)
  6) Downstream assembly metrics + visualization

All major stages are explicitly decoupled:
  - QC once → reuse
  - Assembly without repeating QC
  - smORFs independently
  - AMP prediction independently
  - Refinement independently
  - Downstream analysis independently

The pipeline integrates compositional differential abundance analysis
(ALDEx2) to identify environment-associated smORFs while accounting for
metagenomic count structure and library size variability.
This design enables reproducibility, re-entrancy, and HPC-scalable workflows.
</pre>

<pre>
PIPELINE OVERVIEW
-----------------

Raw reads
   │
   ▼
Quality Control (FastQC / NanoPlot / SeqKit)
   │
   ▼
MetaFlye assembly (--meta)
   │
   ▼
Contig classification (Tiara)
   │
   ├── Bacterial / Prokaryotic contigs
   │      │
   │      ├─ Prodigal → protein-coding ORFs
   │      └─ SmORFinder → smORF candidates
   │
   └── Fungal / Eukaryotic contigs
          │
          └─ MetaEuk → predicted proteins

All predicted peptides
   │
   ▼
peptides_all.faa
   │
   ▼
GLOBAL peptide clustering (MMseqs2)
   │
   ▼
cluster_map.tsv
   │
   ▼
GLOBAL representative CDS catalog
   │
   ▼
cds_nr_global_rep.fna
   │
   ▼
Read mapping per library (minimap2)
   │
   ▼
Primary alignments (MAPQ ≥ 20)
   │
   ▼
Per-CDS abundance tables
   │
   ▼
Compositional normalization (CLR)
   │
   ▼
Differential abundance analysis (ALDEx2)
   │
   ▼
Environment-associated smORFs
   │
   ├─ Bacterial genomic-context refinement
   │
   └─ Fungal genomic-context refinement
        │
        ▼
Refined smORF catalog
   │
   ▼
Macrel AMP prediction
   │
   ▼
predicted_smorfs.with_macrel.tsv
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
   run_smorfs_pipeline.sh                  - smORF + protein discovery pipeline
   smorfs_job.sh                           - Slurm wrapper for smORFs
   metaeuk_job.sh                          - Slurm wrapper for MetaEuk annotation
   run_refine_annot_smorf_bacs.sh          - Bacterial refinement stage
   refine_bacs_job.sh                      - Slurm wrapper for refinement
   refine_bacs.py                          - Genomic-context refinement logic
   downstream_analysis.sh                  - Flye metrics + AMP utilities
   summarize_flye_logs.py                  - Parse Flye logs into metrics table
   plot_metaflye_metrics.R                 - Boxplots from metrics TSV
   predict_amps.py                         - Macrel-based AMP prediction
   attach_macrel_to_predicted_smorfs.py    - Safe Macrel merge (ID-validated)
   run_mmseqs_global_cluster.sh              - Global peptide clustering
   build_global_rep_cds_from_cluster_map.py  - Build NR CDS reference
   run_map_global_cds_abundance.sh          - Runner for GLOBAL CDS mapping step
   map_global_cds_array_task.sh             - Per-library mapping task (minimap2)
   aldex2_global_da.R                     - ALDEx2 differential abundance workflow
   run_differential_abundance_aldex2.sh   - Slurm runner for compositional DA

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
   assembly_metaflye/
   smorfs/
     GLOBAL/
       catalog/
         peptides_all_global.faa
         cds_nr_global_rep.fna
         cds_nr_global_rep.metadata.tsv
       mmseqs/
         cluster_map.tsv
         clusters_pairs.tsv
   abundance/
      global_cds/
         reference/
         manifests/
         aldex2/
            aldex2_counts_matrix.tsv
            aldex2_metadata.tsv
            aldex2_results_kw.tsv
            mapping_summary_plot.pdf
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
 - Global NR CDS catalog ensures consistent abundance units across environments
 - Per-library mapping avoids cross-sample coverage bias
 - Primary alignment filtering (MAPQ ≥ 20)
 - Scratch-aware BAM storage to avoid project directory bloat
 - Manifest auto-generated from FASTQs present in data/
 - Compositional-aware differential abundance using CLR transformation
 - Monte Carlo Dirichlet sampling to model count uncertainty
 - Environment-level statistical comparison across ecosystems
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
  cluster_id refers to the representative ID from MMseqs clustering.
  member IDs are environment-prefixed (ENV|feature_id).

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

<pre>
FUNGAL / EUKARYOTIC smORF REFINEMENT
------------------------------------

Fungal/eukaryotic candidates derived from MetaEuk predictions
can also undergo a refinement stage analogous to bacterial smORF
context evaluation.

This stage integrates:

  • MetaEuk GFF annotations
  • Fungal contig FASTA
  • Global MMseqs recurrence statistics

Three refinement layers are evaluated:

Step 1 — Contig-edge artifact detection
  • Coordinates are parsed from MetaEuk feature_id
  • dist_left / dist_right
  • flag_edge

Step 2 — Overlap with MetaEuk CDS annotations
  • host_cds_id
  • flag_overlap_fraction
  • host_cds_id_embedded
  • flag_embedded

Step 3 — Global structural recurrence
  • cluster_id
  • cluster_size
  • cluster_env_count
  • flag_cluster_recurrent
  • flag_cross_environment

Unlike bacterial refinement, fungal refinement does not rely on
Prodigal annotations and instead uses MetaEuk-derived gene models.

All three steps are independently runnable through the main
workflow launcher (runall.sh).
</pre>

<pre>
EUKARYOTIC / FUNGAL PROTEIN DISCOVERY
-------------------------------------

Eukaryotic contigs are processed separately from bacterial contigs.

Instead of gene modeling pipelines designed for complete genomes
(e.g. Funannotate), the pipeline uses MetaEuk, which is designed
for fragmented metagenomic assemblies.

Workflow:

  fungi_contigs.fasta
       │
       ▼
  MetaEuk prediction
       │
       ├─ metaeuk_preds.fas          (predicted proteins)
       └─ metaeuk_preds.codon.fas    (coding sequences)

Predicted proteins are integrated into:

  peptides_all.faa
  predicted_smorfs.tsv

MetaEuk feature identifiers encode genomic coordinates directly.

Example:

  UniRef50_A0A820GSM0|contig_6867|-|74|1.83e-12|1|4211|4321|...

From this identifier the pipeline extracts:

  contig = contig_6867
  strand = -
  start  = 4211
  end    = 4321

These coordinates are used during fungal refinement to evaluate
contig-edge artifacts and overlaps with MetaEuk CDS annotations.

Because fungal smORF discovery tools are not yet mature,
short proteins are treated as smORF candidates if:

  peptide_length ≤ 100 aa

These candidates are subsequently evaluated by:

  • Macrel AMP prediction
  • global recurrence analysis
  • ecological distribution
</pre>

<pre>
GLOBAL CLUSTERING & RECURRENCE ANALYSIS
----------------------------------------

The pipeline supports global peptide clustering across all environments
to identify structurally recurrent peptide features across ecosystems.

Clustering workflow:

  1) Concatenate peptides_all.faa from:
       RIPARIA_GLOBAL
       PENEIRA_GLOBAL
       CAMPINA_GLOBAL
       UNIFORME_GLOBAL

  2) Prefix FASTA headers with environment label:
       ENV|feature_id

  3) Cluster with MMseqs2:
       --min-seq-id 0.95
       -c 0.8
       --cov-mode 1

  4) Generate:
       clusters_pairs.tsv     (rep <tab> member)
       cluster_map.tsv        (member_id <tab> representative_id)

Refinement Step 3 integrates this information into each smORF catalog:

  cluster_id              Representative sequence ID
  cluster_size            Total members in cluster
  cluster_env_count       Number of distinct environments
  flag_cluster_recurrent  cluster_size ≥ 2
  flag_cross_environment  cluster_env_count ≥ 2

This enables:

  • Detection of globally conserved smORFs
  • Identification of cross-ecosystem recurrence
  • Separation of singletons vs shared features
  • Structured novelty assessment

Clustering can be run independently:

  bash workflow/runall.sh --mmseqs-global-only

Cluster-only refinement mode:

  bash workflow/runall.sh \
    --refine-bacs-only \
    --cluster-only \
    --refine-bacs-sample PENEIRA_GLOBAL

Note: MMseqs clustering is performed on the global pooled peptide catalog
(peptides_all_global.faa). This catalog may contain peptides originating
from both bacterial smORF predictions and fungal/eukaryotic MetaEuk
proteins. Therefore recurrence statistics represent global peptide
similarity patterns rather than fungi-only clustering.
</pre>

<pre>
GLOBAL REPRESENTATIVE CDS CATALOG (NR CDS REFERENCE)
---------------------------------------------------

After global peptide clustering (MMseqs2), the pipeline constructs a
nonredundant nucleotide CDS reference corresponding to the representative
peptide sequences of each cluster.

This step prepares the reference catalog used for downstream abundance
quantification and differential analysis.

Workflow logic:

  1) Read cluster_map.tsv
       (member_id <tab> representative_id)

  2) Extract representative IDs
       ENV|feature_id

  3) Parse ENV prefix to identify source catalog:
       results/smorfs/[ENV]_GLOBAL/catalog/cds_all.fna

  4) Retrieve nucleotide CDS sequence for each representative feature

  5) Build global nonredundant CDS FASTA:

       results/smorfs/GLOBAL/catalog/cds_nr_global_rep.fna

  6) Build accompanying metadata table:

       results/smorfs/GLOBAL/catalog/cds_nr_global_rep.metadata.tsv


FASTA header format:

  ENV|feature_id

Example:

  >UNIFORME|contig_189672_21 # 123 # 456 # 1 # ID=189672_21


Rationale:

Peptide clustering defines functional redundancy.
CDS representatives provide a stable nucleotide reference for
read mapping and abundance estimation.

Clustering is performed only once at the peptide level.
CDS sequences are extracted for the cluster representatives,
ensuring:

  • consistent biological units across environments
  • reduced redundancy
  • compatibility with read-mapping tools
  • reproducible feature identifiers
  • compatibility with DESeq2 count matrices


This step does NOT perform read mapping yet.
It only prepares the global reference catalog.
</pre>

<pre>
READ MAPPING TO GLOBAL CDS CATALOG
----------------------------------

After building the GLOBAL representative CDS catalog, each sequencing
library is mapped independently back to the nonredundant CDS reference.

This produces per-library abundance estimates for each representative smORF/CDS.

Key design principles:

  • mapping performed per library (not per environment)
  • consistent reference across ecosystems
  • primary alignments only
  • MAPQ ≥ 20 filtering
  • scratch-aware BAM storage
  • ALDEx2-compatible compositional count matrices

Mapping workflow:

  data/*.fastq.gz
        │
        ▼
  minimap2 (map-ont)
        │
        ▼
  BAM sorting (temporary)
        │
        ▼
  filtering:
     primary alignments only
     MAPQ ≥ 20
        │
        ▼
  primary_q20.bam
        │
        ▼
  idxstats tables
        │
        ▼
  abundance matrix


Outputs written to scratch:

  /scratch/t.sousa/data_used/read_mapping/

Directory structure:

  read_mapping/
     bam/
        <library_id>/
           *.primary_q20.bam
           *.primary_q20.bam.bai

     stats/
        <library_id>/
           *.primary_q20.flagstat.txt
           *.primary_q20.idxstats.tsv


Representative CDS reference:

  results/abundance/global_cds/reference/global_rep_cds.fna


Manifest automatically built from FASTQs present in:

  data/


Each FASTQ file is treated as one independent library.


Example manifest entry:

  nanopore_shotgun_RDS18_TS-300-mistura_09
      → /scratch/t.sousa/data_used/nanopore_shotgun_RDS18_TS-300-mistura_09.fastq.gz


Rationale:

Mapping reads to the global representative CDS catalog provides
a consistent feature space across all environments.

Counts correspond to structurally nonredundant peptide clusters
defined by MMseqs2.


These abundance estimates support:

  • differential abundance analysis (DESeq2)
  • ecological distribution analysis
  • recurrence validation
  • prioritization of conserved AMP candidates
</pre>

<pre>
DIFFERENTIAL ABUNDANCE ANALYSIS (ALDEx2)
----------------------------------------

After read mapping against the GLOBAL representative CDS catalog,
the pipeline supports compositional differential abundance analysis
using ALDEx2.

ALDEx2 is specifically designed for compositional sequencing data,
addressing biases introduced by library size variation and the
closed-sum constraint inherent to metagenomic count data.

Key characteristics:

  • compositional data aware (CLR transformation)
  • Monte Carlo Dirichlet sampling of count uncertainty
  • robust to uneven sequencing depth
  • suitable for sparse metagenomic feature tables
  • nonparametric statistical testing

Input files:

  /scratch/t.sousa/data_used/read_mapping/stats/<ENV>/<SampleID>/
      *.primary_q20.idxstats.tsv

Each idxstats file provides:

  CDS_id
  CDS_length
  mapped_reads

These values are combined across libraries to produce a global
count matrix with dimensions:

  CDS features × libraries


Statistical workflow:

Step 1 — Feature filtering
  remove CDS with insufficient counts across libraries

Step 2 — Monte Carlo Dirichlet sampling
  models technical variation in observed counts

Step 3 — CLR transformation
  centered log-ratio normalization of compositional data

Step 4 — Global statistical testing
  Kruskal-Wallis test across environments:

    uniforme
    riparia
    peneira
    campina

Step 5 — multiple testing correction
  Benjamini–Hochberg FDR control


Output files:

  results/abundance/global_cds/aldex2/

     aldex2_counts_matrix.tsv
     aldex2_metadata.tsv
     aldex2_results_kw.tsv
     aldex2_significant.tsv
     mapping_summary_plot.pdf


Example output columns:

  feature_id
  kw.ep
  kw.eBH
  effect
  median_uniforme
  median_riparia
  median_peneira
  median_campina


Interpretation:

Significant CDS represent smORFs whose relative abundance
differs across environments after accounting for compositional effects.

These features can be integrated with:

  • recurrence statistics (MMseqs clusters)
  • genomic-context refinement
  • AMP probability scores (Macrel)

enabling prioritization of ecologically relevant antimicrobial candidates.
</pre>

<pre>
REFINEMENT LOGIC OVERVIEW
-------------------------

Refinement (bacterial and fungal) proceeds in structured layers:

Step 1 — Contig-edge artifact detection
  • flag_edge
  • dist_left / dist_right

Step 2 — Embedded / overlap analysis
  • flag_embedded
  • host_cds_id
  • flag_overlap_fraction

Step 3 — Global structural recurrence
  • cluster_id
  • cluster_size
  • cluster_env_count
  • flag_cluster_recurrent
  • flag_cross_environment

These layers allow progressive filtering:

  Structural artifact → Genomic overlap → Ecological recurrence

Each step is independently runnable and does not require re-execution of prior stages.
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

  amp_pred        predicted antimicrobial peptide
  amp_prob        probability of AMP activity
  hemo_pred       predicted hemolytic activity
  hemo_prob       probability of hemolysis
  macrel_class    structural AMP class (if detected)
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

Refinement (Fungal / Eukaryotic):
  --refine-euks-only
  --run-refine-euks
  --refine-euks-sample STR
  --refine-euks-samples FILE
  --refine-euks-step1-only
  --refine-euks-step2-only
  --refine-euks-step3-only
  --refine-euks-steps STR

MetaEuk (Eukaryotic Proteins):
  --metaeuk-only
  --run-metaeuk
  --metaeuk-sample STR
  --metaeuk-samples FILE
  --metaeuk-db PATH

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

GLOBAL CDS ABUNDANCE MAPPING:

Build GLOBAL representative CDS FASTA only:

bash workflow/runall.sh \
  --map-global-cds-build-ref-only


Run read mapping only (libraries currently present in data/):

bash workflow/runall.sh \
  --map-global-cds-only


Build reference + run mapping:

bash workflow/runall.sh \
  --map-global-cds


Restrict mapping to subset of libraries:

bash workflow/runall.sh \
  --map-global-cds-only \
  --map-global-cds-sample-id RDS18


Optional:

--abund-env-name STR
   name of conda env used for mapping tools
   default: smorf_abundance_env

DIFFERENTIAL ABUNDANCE (ALDEx2):

Run full compositional DA workflow:

bash workflow/runall.sh \
  --aldex2-da-only


Run only count matrix preparation:

bash workflow/runall.sh \
  --aldex2-prepare-only


Run only statistical testing:

bash workflow/runall.sh \
  --aldex2-only


Check or install ALDEx2 in abundance environment:

bash workflow/runall.sh \
  --aldex2-check-install-only


Optional parameters:

--aldex2-mc-samples INT
   Monte Carlo instances for Dirichlet sampling
   default: 128

--aldex2-min-count INT
   minimum count threshold per CDS
   default: 10

--aldex2-min-samples INT
   minimum number of libraries passing filter
   default: 2
</pre>

---

<pre>
SCRATCH USAGE
-------------

Large intermediate files are written to scratch to avoid
inflating the project directory.

MetaFlye temporary files:

  /scratch/t.sousa/data_used/metaflye_tmp/<SampleID>/


GLOBAL CDS read mapping outputs:

  /scratch/t.sousa/data_used/read_mapping/

     bam/
        per-library primary alignment BAMs

     stats/
        flagstat and idxstats summaries

     tmp/
        temporary sorted BAM files (auto-deleted)


Final catalogs and metadata remain inside the project directory:

  results/smorfs/
  results/abundance/
</pre>

---

<pre>
EXECUTION EXAMPLES
------------------

QC & ASSEMBLY
-------------

# 1) QC only (default; read-only)
bash workflow/runall.sh

# 2) QC + per-sample MetaFlye assemblies
bash workflow/runall.sh --qc-and-metaflye

# 3) Global co-assembly across all libraries
bash workflow/runall.sh --metaflye-only --global

# 4) Structured global assembly → smORFs → refinement
bash workflow/runall.sh \
  --metaflye-only --global \
  --run-smorfs \
  --run-refine-bacs


smORF DISCOVERY
---------------

# 5) Run smORFs for one assembly
bash workflow/runall.sh \
  --smorfs-only \
  --smorfs-sample TS-0500

# 6) Run smORFs for multiple assemblies
bash workflow/runall.sh \
  --smorfs-only \
  --smorfs-samples metadata/sample_ids.txt


AMP ANNOTATION (Macrel)
-----------------------

# 7) Annotate a predicted smORF catalog with Macrel (ID-safe merge)
bash workflow/runall.sh \
  --macrel-attach-only \
  --predicted-smorfs results/smorfs/TS-0500/catalog/predicted_smorfs.tsv \
  --id-col feature_id \
  --seq-col aa_seq

Macrel fields added:
  amp_pred, amp_prob
  hemo_pred, hemo_prob
  macrel_class


GLOBAL CLUSTERING (Recurrence Analysis)
---------------------------------------

# 8) Build global peptide catalog and cluster across ecosystems
bash workflow/runall.sh \
  --mmseqs-global-only \
  --cpus 16

Outputs:
  results/smorfs/GLOBAL/catalog/peptides_all_global.faa
  results/smorfs/GLOBAL/mmseqs/cluster_map.tsv


BACTERIAL REFINEMENT
--------------------

# 9) Full genomic-context refinement (edge + embedded + recurrence)
bash workflow/runall.sh \
  --refine-bacs-only \
  --refine-bacs-sample PENEIRA_GLOBAL

# 10) Attach recurrence statistics only (Step 3)
bash workflow/runall.sh \
  --refine-bacs-only \
  --cluster-only \
  --refine-bacs-sample PENEIRA_GLOBAL

# 11) Run refinement for multiple ecosystems
bash workflow/runall.sh \
  --refine-bacs-only \
  --refine-bacs-samples metadata/sample_ids.txt


DOWNSTREAM ANALYSIS
-------------------

# 12) Parse Flye logs + generate assembly metrics plots
bash workflow/runall.sh --downstream-only

# 13) Full downstream (metrics + AMP utilities)
bash workflow/runall.sh --downstream-full

FUNGAL REFINEMENT
-----------------

# 14) Full fungal refinement (Steps 1–3)
bash workflow/runall.sh \
  --refine-euks-only \
  --refine-euks-sample PENEIRA_GLOBAL

# 15) Run fungal Step 1 only
bash workflow/runall.sh \
  --refine-euks-only \
  --refine-euks-step1-only \
  --refine-euks-sample PENEIRA_GLOBAL

# 16) Run fungal Step 2 only
bash workflow/runall.sh \
  --refine-euks-only \
  --refine-euks-step2-only \
  --refine-euks-sample PENEIRA_GLOBAL

# 17) Run fungal Step 3 only
bash workflow/runall.sh \
  --refine-euks-only \
  --refine-euks-step3-only \
  --refine-euks-sample PENEIRA_GLOBAL
  
GLOBAL NR CDS REFERENCE
-----------------------

# 18) Build GLOBAL representative CDS catalog from MMseq clusters
bash workflow/runall.sh \
  --build-global-rep-cds

# 19) Run full differential abundance workflow
#     (prepare count matrix + run ALDEx2)

bash workflow/runall.sh \
  --aldex2-da-only

# 20) Check or install ALDEx2 in abundance environment

bash workflow/runall.sh \
  --aldex2-check-install-only


# 21) Generate mapping summary plots only
#     (mapped %, primary mapped %, MAPQ ≥ 20 filtering)

bash workflow/runall.sh \
  --aldex2-flagstat-only


# 22) Prepare count matrix + metadata only

bash workflow/runall.sh \
  --aldex2-prepare-only


# 23) Run statistical testing only
#     (requires prepared count matrix)

bash workflow/runall.sh \
  --aldex2-only

# 24) Increase Monte Carlo precision

bash workflow/runall.sh \
  --aldex2-da-only \
  --aldex2-mc-samples 256

# 25) Modify filtering thresholds

bash workflow/runall.sh \
  --aldex2-da-only \
  --aldex2-min-count 20 \
  --aldex2-min-samples 3

# 26) Specify environments explicitly

bash workflow/runall.sh \
  --aldex2-da-only \
  --aldex2-envs uniforme,riparia,peneira,campina

# 27) Custom output directory

bash workflow/runall.sh \
  --aldex2-da-only \
  --aldex2-outdir results/custom_aldex2
</pre>

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
