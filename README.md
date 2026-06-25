# MGall_Analysis

Comprehensive bioinformatics pipeline for *Mycoplasma gallisepticum* genomic analysis, including genome assembly, annotation, phylogenetics, pangenome analysis, GWAS, epidemiological modeling, and advanced phylodynamic inference.

**Author**: Dr. Alberto Torcello-Requena  
**Affiliation**: Environment and Sustainability Institute, University of Exeter  
**Contact**: <a.torcello-requena@exeter.ac.uk>

---

## Table of Contents

1. Aim & Description
2. Repository Structure
3. Installation & Requirements
4. Workflow Overview
5. Detailed Script Documentation
6. R Analysis Pipeline
7. Typical Use Cases
8. Troubleshooting
9. Citation

---

## Aim & Description

### Project Objectives

This repository provides an end-to-end computational framework for comprehensive genomic analysis of *Mycoplasma gallisepticum* (*Mg*) sampled from wild house finches (*Haemorhous mexicanus*). The pipeline integrates:

- **Population-level inference**: tracking lineage emergence, fitness trajectories, and population dynamics
- **Genomic characterization**: identifying lineage-defining mutations, virulence-associated variants, and adaptive signatures
- **Synteny & virulence context**: mapping GWAS-significant loci within genomic neighborhoods to infer functional impacts
- **Epidemiological prediction**: parameterizing temporal population models from phylogenetic data

### Scientific Context

*Mycoplasma gallisepticum* (*Mg*) is a pathogenic bacterium causing conjunctivitis and respiratory disease in wild and domestic birds. Understanding its genomic evolution, population structure, and virulence determinants is critical for:

- Predicting outbreak dynamics in wild bird populations
- Identifying lineages with increased transmission or virulence
- Uncovering molecular mechanisms of disease severity
- Guiding surveillance and intervention strategies

This pipeline supports research on host-pathogen interactions, evolutionary fitness, and disease ecology by integrating phylogenetic, genomic, and epidemiological data.

---

## Overview

This repository contains a complete computational pipeline for analyzing *Mycoplasma gallisepticum* genomic data from house finches (*Haemorhous mexicanus*). The analysis integrates:

- **Data management**: Sample manifest generation, metadata extraction, MD5 validation, SRA data retrieval
- **Genome assembly**: De novo assembly (SPAdes), quality control with read trimming (Sickle)
- **Annotation**: Prokka-based structural annotation, protein extraction, gene context analysis
- **Alignment & filtering**: Sequence alignment trimming, gap removal, parsimony-informative site selection
- **Phylogenetics**: Maximum likelihood trees (IQ-TREE), time-calibrated phylogenies (BEAST), recombination detection (Gubbins), topology hypothesis testing (SH/AU tests)
- **Pangenomics**: Gene presence/absence analysis (Panaroo), cluster characterization, lineage-specific gene identification, synteny visualization
- **Population genetics**: SNP calling (bcftools), GWAS analysis, Manhattan plots with gene annotation, dN/dS analysis
- **Virulence analysis**: Phage detection (VIBRANT), defense system annotation (PADLOC), motility gene characterization
- **Phylodynamics**: Bayesian skyline/skygrid plots, lineage fitness estimation, phylogenetic signal analysis, Brownian motion trait evolution
- **Epidemiological modeling**: Temporal dynamics, effective population size inference, phylowave analysis for selective sweeps

---

## Repository Structure

```bash
MGall_Analysis/
├── .git/                           # Version control
├── .gitignore                      # Git ignore rules
├── .idea/                          # PyCharm project settings
├── .vscode/                        # VS Code settings
├── README.md                       # This file
├── MGall_R_analysis/               # R-based analyses & phylodynamics
│   ├── 2_Functions/                # Reusable R utility functions
│   │   ├── 2_1_Index_computation_20240909.R
│   │   ├── 2_2_Lineage_detection_20240909.R
│   │   ├── 2_3_Lineage_fitness_20240909.R
│   │   └── 2_4_Lineage_defining_mutations.R
│   ├── Brownian_Alberto/           # Brownian motion trait evolution
│   │   ├── Brownian_Alberto.R      # Main analysis script
│   │   ├── Brownian_Alberto_original.R
│   │   ├── Brownian_Alberto_Lineages.R
│   │   ├── coretree_noR.nwk        # Phylogenetic tree (Newick)
│   │   ├── Edited_VA94_consensus_all_trimmed_60threshold_50_combined.finaltree.newick
│   │   ├── leaves_L1.txt           # Lineage 1 sample IDs
│   │   ├── leaves_L2.txt           # Lineage 2 sample IDs
│   │   ├── meta_original.txt       # Original metadata
│   │   └── meta.txt                # Sample metadata (dates, locations, traits)
│   ├── BacDating.R                 # Bayesian dating analysis (divergence times)
│   ├── Metadata_trimming.py        # Python utility for metadata processing
│   ├── MGall_R_analysis.Rproj      # RStudio project file
│   ├── Panstripe.R                 # Pangenome stripe plots with bootstrap CIs
│   ├── Phylowave.R                 # Phylogenetic wave analysis (selective sweeps)
│   ├── Position_segregating_sites.R # Segregating site analysis
│   ├── Skyline_lineages.R          # Skyline plots per lineage
│   ├── Skylineplot.R               # BEAST skyline plot visualization
│   └── SNPs_tree_making.r          # SNP matrix tree construction
├── Scripts/                        # Python & shell scripts
│   ├── General_python_scripts/     # General-purpose utilities
│   │   ├── Alignment2Fasta.py      # Convert VCF/other formats to FASTA
│   │   ├── Alignment_trimming.py   # Remove gap-only columns
│   │   ├── Alignment_trimming_informativesites.py  # Keep parsimony-informative sites only
│   │   ├── Alignment_trimming_threshold.py  # Filter by sequence identity threshold
│   │   ├── AT_ratio_overtime.py    # AT content over time regression
│   │   ├── GFFs2FNA.py             # Extract sequences from GFF annotations
│   │   ├── GeneContext_verification.py  # Validate gene neighborhood integrity
│   │   ├── Pangenome_translation.py  # Translate protein sequences
│   │   ├── Position_segregating_sites.py  # Extract segregating sites from BAM
│   │   ├── Preparation_trees.py    # Prepare tree files for analysis
│   │   ├── Rename_dictionary.py    # Create name replacement dictionaries
│   │   ├── Rename_files.py         # Batch file renaming
│   │   ├── Rename_GFFs_metadata.py # Rename GFF coordinate systems
│   │   ├── Rename_taxa_columns.py  # Rename dataframe columns
│   │   ├── Rename_taxa_file.py     # Rename sequences in files
│   │   ├── Rename_taxa_folders.py  # Rename directory structures
│   │   └── Tree_VCF.py             # Build alignment from VCF files
│   ├── MG_genome/                  # Core genomic analyses
│   │   ├── Alignment_gubbins_comparison.py  # Compare pre/post-Gubbins alignments
│   │   ├── Alignment_lineages.py   # Lineage-specific alignment extraction
│   │   ├── Analysis_BEAST_bifurcation.py  # Bifurcation timing from BEAST trees
│   │   ├── Ancestral_tree_metadata.py  # Annotate tree with ancestral states
│   │   ├── Check_GM_context_COG_mortality.py  # Synteny context (virulence phenotype)
│   │   ├── Check_GM_context_COG_swelling.py   # Synteny context (alternative phenotype)
│   │   ├── Creation_list_input_gubbins.py  # Prepare file lists for Gubbins
│   │   ├── Dictionary_clusters_heatmap.py  # Categorize pangenome clusters (virulence/motility/etc.)
│   │   ├── dnds_data_preparation.py  # Prepare dN/dS calculation data
│   │   ├── GWAS_results_fasta_results.py  # Extract GWAS-significant gene sequences
│   │   ├── Lineage_differences.py  # Identify lineage-specific genes
│   │   ├── Lineage_differences_vcf.py  # SNP-level lineage differences
│   │   ├── Manhattan_annotated.py  # Publication-quality Manhattan plots
│   │   ├── Manhattan_plot.py       # Basic Manhattan plot generation
│   │   ├── Manhattan_prep.py       # Prepare GWAS data for plotting
│   │   ├── Manifest_files.py       # Generate NCBI submission manifest
│   │   ├── Metadata_extraction.py  # Extract metadata from GenBank/SRA
│   │   ├── Metadata_Viral_analysis.py  # Phage/prophage metadata integration
│   │   ├── Pangenome_clustering_heatmap.py  # Hierarchical clustering visualization
│   │   ├── Pangenome_extraction.py  # Extract genes by presence/absence pattern
│   │   ├── Pastml_analysis.py      # Ancestral state reconstruction (geographic/phenotypic)
│   │   ├── Phylogeny_snps.py       # Add SNP distribution nodes to trees
│   │   ├── Plot_GM_context.py      # Plot synteny/context with GWAS annotation
│   │   ├── Plot_metadata.py        # Metadata distribution plots
│   │   ├── Presence_absence_matrix.py  # Build gene presence/absence table
│   │   ├── SkyGrid_plot.py         # BEAST Skygrid effective population size plots
│   │   ├── SkyGrid_tree.py         # Skygrid tree preparation
│   │   └── Tree_lineages.py        # Assign lineage IDs to tree tips
│   ├── MG_motility/                # Motility gene-specific analysis
│   │   ├── Analysis_vcf_motility.py  # SNP analysis in motility-related genes
│   │   └── Motility_protein_preparation.py  # Motility protein database creation
│   ├── Shell_scripts/              # Bash pipeline scripts for HPC/local execution
│   │   ├── BCFtools.sh             # Consensus sequence generation from BAM
│   │   ├── Download_sra.sh         # Batch download FASTQ from EBI/SRA
│   │   ├── Gubbins_pipeline.sh     # Recombination detection and masking
│   │   ├── Lineages_states.sh      # Lineage state inference pipeline
│   │   ├── Lucy_Prokka.sh          # Prokka annotation wrapper
│   │   ├── Lucy_SPADES_Assembly.sh # SPAdes assembly wrapper
│   │   ├── Padloc.sh               # Defense system annotation (CRISPR/etc.)
│   │   ├── Panaroo_pipeline.sh     # Pangenome construction
│   │   ├── Poppunk_code.sh         # PopPUNK clustering (alternative lineage detection)
│   │   ├── Prokka_command.sh       # Individual Prokka annotation runs
│   │   ├── Prokka_complete_genomes.sh  # Prokka for complete/finished genomes
│   │   ├── SH_comparison.sh        # Tree topology hypothesis tests
│   │   ├── Sickle_readTrimming.sh  # Quality trimming (paired-end FASTQ)
│   │   ├── SPADES_Pipeline.sh      # SPAdes assembly pipeline
│   │   ├── Tree_formation.sh       # Phylogenetic tree file preparation
│   │   ├── Treemaking_commands.sh  # IQ-TREE command generation and execution
│   │   └── Vibrant.sh              # Prophage/phage prediction
│   └── Slurm_scripts/              # HPC cluster job submission (SLURM)
│       ├── IQ-Tree_slurm.sh        # ML phylogeny inference (scalable)
│       ├── Mgall_mapping_slurm.sh  # Read mapping to reference (BBMap)
│       └── Variant_calling.sh      # SNP/indel calling (bcftools)
```

---

## Installation & Requirements

### Python Environment

**Recommended**: Python 3.8+

```bash
# Core packages
pip install numpy pandas scipy matplotlib seaborn biopython bcbio-gff openpyxl tqdm

# Optional: Interactive plotting
pip install PyQt5  # or PySide6
pip install ipympl  # for Jupyter

# Conda alternative
conda create -n MGall_Analysis python=3.10 \
  numpy pandas scipy matplotlib seaborn biopython openpyxl tqdm -c conda-forge
conda activate MGall_Analysis
```

### R Environment

**Required R packages**:

```r
install.packages(c("ape", "phangorn", "BactDating", "phytools", "coda", 
                   "vcfR", "ggplot2", "ggtree", "cowplot", "mgcv", 
                   "devtools", "RColorBrewer"))

devtools::install_github("laduplessis/bdskytools")
devtools::install_github("laduplessis/beastio")
```

### External Tools

Required command-line tools (install via conda/apt/homebrew):

- **Assembly**: SPAdes, sickle
- **Mapping**: BBMap/BBTools, samtools, bcftools
- **Annotation**: Prokka, VIBRANT, PADLOC
- **Phylogenetics**: IQ-TREE, BEAST2, Gubbins, RAxML, FastTree
- **Pangenomics**: Panaroo, PopPUNK
- **Utilities**: seqkit, muscle

---

## Workflow Overview

### 1. Data Preparation

```bash
# Download SRA data
bash Scripts/Shell_scripts/Download_sra.sh

# Quality trimming
bash Scripts/Shell_scripts/Sickle_readTrimming.sh

# Generate manifest
python Scripts/MG_genome/Manifest_files.py
```

### 2. Genome Assembly & Annotation

```bash
# De novo assembly
bash Scripts/Shell_scripts/Lucy_SPADES_Assembly.sh /path/to/reads/

# Structural annotation
bash Scripts/Shell_scripts/Lucy_Prokka.sh /path/to/assemblies/

# Extract genome sequences from GFF
python Scripts/General_python_scripts/GFFs2FNA.py
```

### 3. Mapping & Variant Calling

```bash
# Map to reference (HPC)
sbatch Scripts/Slurm_scripts/Mgall_mapping_slurm.sh

# Call variants
sbatch Scripts/Slurm_scripts/Variant_calling.sh

# Generate consensus sequences
bash Scripts/Shell_scripts/BCFtools.sh
```

### 4. Phylogenetic Analysis

```bash
# Alignment trimming
python Scripts/General_python_scripts/Alignment_trimming.py \
  --sequences alignment.fasta --output trimmed.fasta

# Recombination filtering
bash Scripts/Shell_scripts/Gubbins_pipeline.sh

# Maximum likelihood tree
sbatch Scripts/Slurm_scripts/IQ-Tree_slurm.sh

# Tree topology comparison
bash Scripts/Shell_scripts/SH_comparison.sh
```

### 5. Pangenome Analysis

```bash
# Build pangenome
bash Scripts/Shell_scripts/Panaroo_pipeline.sh

# Cluster analysis
python Scripts/MG_genome/Dictionary_clusters_heatmap.py

# Lineage-specific genes
python Scripts/MG_genome/Lineage_differences.py
```

### 6. GWAS & Association Studies

```bash
# Prepare Manhattan plot data
python Scripts/MG_genome/Manhattan_prep.py

# Generate annotated plot
python Scripts/MG_genome/Manhattan_annotated.py \
  --manhattan gwas.tsv --gff reference.gff --output plot.png

# Extract significant genes
python Scripts/MG_genome/GWAS_results_fasta_results.py
```

#### Synteny/context analysis (Panaroo clusters + GWAS direction)

Scripts:

- Scripts/MG_genome/Check_GM_context.py — builds per-cluster synteny context from Panaroo CSV/Rtab/GML and sample GFFs.
- Scripts/MG_genome/Plot_GM_context.py — plots present and absent contexts, annotates GWAS beta direction, and summarizes signatures.

Outputs per significant cluster:

- cluster_synteny.tsv — per-sample rows with target/product, neighbor products/clusters, presence flags, GWAS beta, direction match, absence classification.
- cluster_signatures.tsv — product signature frequencies.
- cluster_cluster_signatures.tsv — cluster signature frequencies.
- cluster_consensus_context.tsv — consensus neighbor clusters among present samples.
- Plots/PNG: synteny tiles (products/clusters), signature bars, absence summary, consensus neighbors, plus a sample status TSV.

GWAS beta interpretation:

- beta ≥ 0: presence → higher virulence.
- beta < 0: absence → higher virulence.
Tile plots mark center tiles: green=matches direction, red=mismatch, gray=unknown. Absent rows get borders by class (annotation_issue, region_missing, partial_region_loss, target_only_missing_or_misannotated).

Usage:

```bash
# Generate context tables (paths configured inside script)
python Scripts/MG_genome/Check_GM_context.py

# Plot all clusters, include absent rows, enlarge tiles
python Scripts/MG_genome/Plot_GM_context.py \
  --ctx-dir synteny_context_mortality_COGs \
  --include-absent \
  --mode both \
  --tile-col-w 1.3 --tile-row-h 0.45 --dpi 260
```

Troubleshooting (plotting):

- Headless backend: set matplotlib to Agg to avoid Qt/Wayland errors.

```python
import matplotlib
matplotlib.use('Agg')
```

- If figures clip or overlap, increase tile sizes via flags: --tile-col-w, --tile-row-h, and --dpi. Legend space is auto-reserved.

### 7. Temporal & Phylodynamic Analysis

```r
# In R
source("MGall_R_analysis/Skylineplot.R")
source("MGall_R_analysis/BacDating.R")
source("MGall_R_analysis/Phylowave.R")
```

---

## Detailed Script Documentation

### General Python Scripts

#### **Alignment_trimming.py**

Removes gap-only columns from multiple sequence alignments.

```bash
python Alignment_trimming.py --sequences input.fasta --output nogaps.fasta
```

**Outputs**: Trimmed alignment + CSV of removed columns.

#### **Alignment_trimming_informativesites.py**

Filters alignment to parsimony-informative sites only.

```bash
python Alignment_trimming_informativesites.py
```

**Use case**: Reducing alignment size for faster phylogenetic inference.

#### **AT_ratio_overtime.py**

Calculates AT content over time and performs regression analysis.

**Outputs**:

- `at_ratio_vs_time.png`
- OLS regression summary (console)

#### **Position_segregating_sites.py**

Identifies segregating sites and extracts read-level information from BAM files.

**Outputs**: `Test_pos.csv`

#### **Rename_dictionary.py**

Creates standardized name replacement dictionaries for Lucy/SRA samples.

**Outputs**:

- `Lucy_replacements.pickle`
- `Camille_replacements_foldername.pickle`

#### **Rename_taxa_columns.py** / **Rename_taxa_file.py** / **Rename_taxa_folders.py**

Batch renaming utilities for various data types.

```bash
python Rename_taxa_columns.py --file input.tsv        # Rename dataframe columns
python Rename_taxa_file.py --file sequences.fasta     # Rename sequences in FASTA
python Rename_taxa_folders.py --directory /path/to/data/  # Rename directory structures
```

#### **Alignment_trimming_threshold.py**

Filters sequence alignment based on sequence identity threshold.

**Use case**: Remove divergent sequences or enforce minimum sequence similarity.

**Outputs**: Filtered alignment FASTA file.

---

### MG_genome Scripts

#### **Manifest_files.py**

Generates NCBI submission manifest with MD5 checksums.

**Features**:

- Reads `.md5` files or Excel sheet (`md5_values.xlsx`)
- Normalizes FASTQ filename suffixes
- Handles paired-end reads

**Outputs**:

- `MGall_Manifest.csv`
- `MGall_Manifest_with_MD5.csv`

#### **Metadata_extraction.py**

Extracts metadata from GenBank files and combines with SRA/sample data.

**Outputs**: `Metadata_genomes.xlsx` (multi-sheet Excel workbook)

#### **Presence_absence_matrix.py**

Builds gene presence/absence matrix from Panaroo outputs.

**Inputs**: Panaroo CSV and Rtab files

**Outputs**: Binary presence/absence TSV matrix (samples × genes)

#### **Pangenome_extraction.py**

Extracts specific gene clusters based on presence/absence patterns.

**Features**: Filter by cluster ID, presence threshold, lineage-specific genes

**Outputs**: FASTA files of extracted genes

#### **Alignment_gubbins_comparison.py**

Compares sequence alignments before and after Gubbins recombination masking.

**Metrics**: SNP differences, recombined regions, sequence divergence

**Outputs**: Comparison statistics and plots

#### **Alignment_lineages.py**

Extracts lineage-specific alignment subsets from full alignment.

**Workflow**: 
1. Load sample-to-lineage mapping
2. Filter alignment to lineage samples
3. Export lineage-specific FASTA

**Outputs**: Per-lineage alignment FASTA files

#### **Analysis_BEAST_bifurcation.py**

Analyzes bifurcation times from BEAST time-calibrated trees.

**Extracts**: Node ages, lineage divergence times, HPD intervals

**Outputs**: Bifurcation timing table, ancestor-descendant lineage pairs

#### **Ancestral_tree_metadata.py**

Annotates internal tree nodes with ancestral states and metadata.

**Features**: Links internal nodes to reconstructed sequences/traits

**Outputs**: Annotated Newick tree with node metadata

#### **Lineage_differences_vcf.py**

Identifies SNPs that distinguish lineages using VCF data.

**Workflow**:
1. Load VCF files (alignment or individual sample VCFs)
2. Define lineage membership from phylogenetic tree
3. Extract lineage-diagnostic SNPs
4. Classify by fixation status

**Outputs**: SNP matrix, variant classification table

#### **Pastml_analysis.py**

Performs ancestral state reconstruction (ASR) for geographic locations or phenotypic traits.

**Methods**: Maximum likelihood, joint reconstruction

**Features**: Maps trait evolution onto phylogeny

**Outputs**: Annotated tree with ancestral states at internal nodes

#### **Plot_metadata.py**

Generates distribution plots of sample metadata (collection dates, locations, traits).

**Visualization types**: Time series, geographic distribution, trait distributions

**Outputs**: PNG/PDF figures with publication-ready formatting

#### **SkyGrid_tree.py**

Prepares phylogenetic trees for BEAST Skygrid analysis.

**Workflow**: 
1. Validates tree format (Newick)
2. Extracts sample dates/ages
3. Formats for BEAST input

**Outputs**: BEAST-compatible tree file with date metadata

#### **dnds_data_preparation.py**

Prepares codon alignment data for dN/dS ratio calculation.

**Workflow**:
1. Extracts coding sequences from GFF + genome
2. Creates codon alignments
3. Formats for dN/dS tools (e.g., codeml)

**Outputs**: Codon alignment FASTA, control files

#### **Lineage_differences.py**

Identifies lineage-specific genes using presence/absence matrix.

**Key parameters**:

- `filter1`: Minimum presence in lineage 1 (default: 0.60)
- `filter2`: Maximum presence in lineage 2 (default: 0.30)

**Outputs**:

- Venn diagrams
- Heatmaps
- FASTA files of lineage-specific genes

#### **Manhattan_annotated.py**

Creates publication-quality Manhattan plots with gene annotations.

```bash
python Manhattan_annotated.py \
  --manhattan mortality_manhattan.txt \
  --gff reference.gff \
  --output manhattan.png \
  --distance 5000
```

**Features**:

- Bonferroni correction threshold
- Configurable annotation distance window
- Handles GFF coordinate systems

#### **Dictionary_clusters_heatmap.py**

Categorizes pangenome clusters (lipoproteins, virulence, Cas9, motility).

**Outputs**: `cluster_dict.pickle`

#### **Pangenome_clustering_heatmap.py**

Generates hierarchical clustering heatmap of gene presence/absence.

**Features**:

- Lineage-based sample ordering
- Customizable color schemes
- Exports sorted gene lists

#### **SkyGrid_plot.py**

Visualizes BEAST Skygrid output with lineage-specific trajectories.

**Outputs**:

- `Combined_Figures.png` (effective population size + Re)
- Individual plots for debugging

#### **Phylogeny_snps.py**

Adds pie chart nodes to phylogenetic trees showing SNP distributions.

**Outputs**: Annotated tree visualization

---

### MG_motility Scripts

#### **Analysis_vcf_motility.py**

Analyzes SNPs in motility-related genes from VCF files.

**Workflow**:

1. Loads candidate motility genes from BLAST results
2. Filters pangenome for non-motility genes
3. Extracts protein sequences
4. Builds SNP matrix across samples

**Outputs**: `Matrix_mutations_SNPs.tsv`

#### **Motility_protein_preparation.py**

Creates protein database from Prokka annotations for motility analysis.

**Outputs**: `Lucy_protein_db.fna`

---

### Shell Scripts

#### **Download_sra.sh**

Batch download of FASTQ files from EBI FTP.

```bash
bash Download_sra.sh
```

**Features**: Parallel downloads using background processes.

#### **Sickle_readTrimming.sh**

Quality trimming with Sickle (paired-end mode).

```bash
bash Sickle_readTrimming.sh
```

**Outputs**: `*_sickle_R1.fastq.gz`, `*_sickle_R2.fastq.gz`, `*_sickle_single.fastq.gz`

#### **Lucy_SPADES_Assembly.sh**

De novo assembly pipeline for Lucy samples.

```bash
bash Lucy_SPADES_Assembly.sh /path/to/reads/
```

**Skips existing assemblies** to enable resumable runs.

#### **Lucy_Prokka.sh** / **Prokka_command.sh**

Structural annotation with Prokka (genetic code 4 for Mycoplasma).

```bash
bash Lucy_Prokka.sh /path/to/assemblies/
```

**Parameters**:

- Genus: *Mycoplasma*
- Species: *gallisepticum*
- Code: 4 (Mycoplasma genetic code)

#### **Gubbins_pipeline.sh**

Recombination detection and masking.

```bash
bash Gubbins_pipeline.sh
```

**Key options**:

- `--filter-percentage 50`: Remove sites with >50% recombination
- `--extensive-search`: Improved accuracy
- `--recon-with-dates`: Time-aware reconstruction

#### **Panaroo_pipeline.sh**

Pangenome construction with strict cleaning mode.

```bash
bash Panaroo_pipeline.sh
```

**Parameters**:

- `--clean-mode strict`
- `--core_threshold 0.98`
- `--codon-table 4`

#### **SH_comparison.sh**

Tree topology hypothesis testing (SH test, AU test).

```bash
bash SH_comparison.sh
```

**Combines trees from**:

- 60% threshold alignment
- 80% threshold alignment
- Parsimony-informative sites only
- No-gaps alignment

#### **Vibrant.sh** / **Padloc.sh**

Prophage detection and defense system annotation.

```bash
bash Vibrant.sh /path/to/genomes/
bash Padloc.sh /path/to/genomes/
```

#### **SPADES_Pipeline.sh**

De novo assembly pipeline using SPAdes (alternative to Lucy_SPADES_Assembly.sh).

**Parameters**: Auto-detection of read length and insert size

#### **Poppunk_code.sh**

PopPUNK clustering for lineage/strain assignment (alternative phylogenetic approach).

**Features**: Rapid lineage identification, strain-level resolution

**Outputs**: Cluster assignments, pairwise distances

#### **Tree_formation.sh**

Prepares tree files for downstream analyses (format conversion, validation).

**Supports**: Newick, Nexus, PhyloXML formats

**Outputs**: Standardized tree files

#### **Treemaking_commands.sh**

Generates and executes IQ-TREE commands for phylogenetic inference.

**Includes**: Model selection, bootstrap replication, alternative topologies

**Outputs**: Inferred trees with support values

#### **Lineages_states.sh**

Lineage state inference pipeline (trait/phenotype assignment to lineages).

**Workflow**: 
1. Maps samples to lineages using phylogenetic tree
2. Aggregates phenotypic/epidemiological data per lineage
3. Infers lineage characteristics

**Outputs**: Lineage state summary tables

---

### Slurm Scripts (HPC)

#### **Mgall_mapping_slurm.sh**

Maps reads to reference genome using BBMap.

```bash
sbatch Mgall_mapping_slurm.sh
```

**Resource requirements**: 20 CPUs, 23 GB RAM per node.

#### **Variant_calling.sh**

Calls SNPs with bcftools (mpileup + call).

```bash
sbatch Variant_calling.sh
```

**Filters**: Minimum quality score ≥10, indels excluded.

#### **IQ-Tree_slurm.sh**

Maximum likelihood tree inference with model testing.

```bash
sbatch IQ-Tree_slurm.sh
```

**Features**: Automatic model selection, 1000 bootstraps.

---

## R Analysis Pipeline

### Core Functions (2_Functions/)

#### **2_1_Index_computation_20240909.R**

Computes phylogenetic indices and summary statistics.

**Computes**: Tree balance, lineage diversity, node support metrics

**Key outputs**: Index matrix for downstream analyses

#### **2_2_Lineage_detection_20240909.R**

Tools for identifying lineages using phylogenetic clustering.

**Key functions**:

- `number.descendants.all.nodes()`: Compute descendant counts
- `test_node_gam()`: GAM-based lineage significance testing

#### **2_3_Lineage_fitness_20240909.R**

Fitness estimation for detected lineages.

**Outputs**: MCMC chains for fitness parameters.

#### **2_4_Lineage_defining_mutations.R**

Ancestral state reconstruction for SNPs/amino acids.

**Workflow**:

1. Load VCF files for coding regions (E, M, N, ORF1a, etc.)
2. Reconstruct ancestral states using maximum parsimony
3. Identify lineage-defining mutations

**Outputs**: RDS files with reconstruction objects.

---

### Standalone Scripts

#### **BacDating.R**

Bayesian dating of phylogenies using [`BactDating`](https://github.com/xavierdidelot/BactDating) package.

**Use case**: Estimate divergence times and substitution rates.

#### **Skylineplot.R**

Visualizes BEAST Bayesian Skyline/Skygrid analyses.

**Features**:

- HPD interval plotting
- Multiple trajectory comparison
- Time-calibrated Ne(τ) curves

#### **Panstripe.R**

Generates pangenome stripe plots showing gene gain/loss patterns.

```r
source("MGall_R_analysis/Panstripe.R")
```

**Outputs**: PNG files with bootstrap confidence intervals.

#### **Phylowave.R**

Phylogenetic wave analysis for detecting selective sweeps.

**Inputs**: Time-calibrated tree + SNP matrix.

#### **Brownian_Alberto.R**

Tests virulence traits for Brownian motion evolution.

**Features**:

- Phylogenetic signal (λ) estimation
- Ancestral state reconstruction
- Virulence index analysis

**Metadata**: `meta.txt` contains sample dates, locations, virulence scores.

#### **SNPs_tree_making.r**

Constructs phylogenetic trees directly from SNP matrices.

**Workflow**:
1. Loads SNP matrix (samples × variable sites)
2. Computes pairwise distances
3. Builds neighbor-joining or UPGMA tree
4. Visualizes with bootstrap support

**Outputs**: Newick tree file, tree plot

#### **Position_segregating_sites.R**

Analyzes segregating sites across samples in phylogenetic context.

**Features**: 
- Counts shared/unique polymorphisms
- Computes nucleotide diversity
- Maps segregating sites to genes

**Outputs**: Position tables, diversity plots

#### **Skyline_lineages.R**

Generates Bayesian Skyline plots separately for each detected lineage.

**Features**:
- Lineage-specific population size trajectories
- HPD interval shading
- Comparative trajectory visualization

**Outputs**: Per-lineage PNG plots, combined figures



## Typical Use Cases

### Use Case 1: New Sample Processing

```bash
# 1. Download and trim
bash Scripts/Shell_scripts/Download_sra.sh
bash Scripts/Shell_scripts/Sickle_readTrimming.sh

# 2. Assembly and annotation
bash Scripts/Shell_scripts/Lucy_SPADES_Assembly.sh /path/to/reads/
bash Scripts/Shell_scripts/Lucy_Prokka.sh /path/to/assemblies/

# 3. Extract metadata
python Scripts/MG_genome/Metadata_extraction.py

# 4. Generate manifest
python Scripts/MG_genome/Manifest_files.py
```

### Use Case 2: Phylogenetic Analysis

```bash
# 1. Map to reference
sbatch Scripts/Slurm_scripts/Mgall_mapping_slurm.sh

# 2. Call variants
sbatch Scripts/Slurm_scripts/Variant_calling.sh

# 3. Build alignment
python Scripts/General_python_scripts/Tree_VCF.py \
  --directory /path/to/vcf/ --output alignment.fasta

# 4. Trim and filter
python Scripts/General_python_scripts/Alignment_trimming.py \
  --sequences alignment.fasta --output trimmed.fasta

# 5. Detect recombination
bash Scripts/Shell_scripts/Gubbins_pipeline.sh

# 6. Build tree
sbatch Scripts/Slurm_scripts/IQ-Tree_slurm.sh
```

### Use Case 3: GWAS Workflow

```bash
# 1. Prepare association test output
python Scripts/MG_genome/Manhattan_prep.py

# 2. Annotate significant SNPs
python Scripts/MG_genome/Manhattan_annotated.py \
  --manhattan mortality_manhattan.txt \
  --gff VA94_reference.gff \
  --output mortality_manhattan.png

# 3. Extract gene sequences
python Scripts/MG_genome/GWAS_results_fasta_results.py

# 4. Analyze motility genes (if relevant)
python Scripts/MG_motility/Analysis_vcf_motility.py
```

### Use Case 4: Pangenome Analysis

```bash
# 1. Build pangenome
bash Scripts/Shell_scripts/Panaroo_pipeline.sh

# 2. Define gene clusters
python Scripts/MG_genome/Dictionary_clusters_heatmap.py

# 3. Identify lineage-specific genes
python Scripts/MG_genome/Lineage_differences.py

# 4. Visualize in R
Rscript MGall_R_analysis/Panstripe.R
```

### Use Case 5: Temporal Dynamics

```bash
# 1. Run BEAST analysis (external)
# 2. Process Skygrid output
python Scripts/MG_genome/SkyGrid_plot.py

# 3. Visualize in R
Rscript MGall_R_analysis/Skylineplot.R
```

---

## Troubleshooting

### Common Issues

#### **1. MD5 Checksum Mismatch**

```bash
Error: MD5 file contains extra text
```

**Solution**: Ensure `.md5` files contain only the 32-character hash:

```bash
md5sum file.fastq > file.fastq.md5
```

#### **2. GFF Annotation Failures**

```bash
Unknown gene (0bp away)
```

**Solutions**:

- Check chromosome name matching: `grep "^>" reference.fna` vs `awk '{print $1}' reference.gff`
- Verify coordinate system (GFF3 is 1-based)
- Increase distance window: `--distance 10000`

#### **3. Biopython Deprecation Warnings**

```python
BiopythonDeprecationWarning: feature.strand is deprecated
```

**Solution**: Use `feature.location.strand` instead.

#### **4. Qt Backend Errors**

```python
qt.qpa.plugin: Could not load the Qt platform plugin "xcb"
```

**Solution**: Set matplotlib backend in script header:

```python
import matplotlib
matplotlib.use('Agg')
```

#### **5. Prokka Fails on Large Genomes**

```bash
Error: Too many contigs
```

**Solution**: Filter short contigs before annotation:

```bash
seqkit seq -m 500 assembly.fasta > filtered.fasta
```

#### **6. IQ-TREE Memory Issues**

```bash
ERROR: Not enough memory
```

**Solution**: Reduce alignment size or use simpler model:

```bash
iqtree -s alignment.fasta -m GTR+G -mem 10G
```

---

### Data File Locations

Key hardcoded paths in scripts (update for your system):

```bash
# Base directory
/home/albertotr/OneDrive/Data/Cambridge_Project/

# Subdirectories
Lucy_reads/            # Raw Lucy sequencing data
SRA_reads/             # Downloaded SRA data
CheckMbins/            # Reference genomes
Mapped_output_*/       # Mapping results
pangenome_results_*/   # Panaroo outputs
GWAS/                  # Association study results
```

**Recommendation**: Use environment variables or config files for production use.

---

## Citation

If you use this pipeline, please cite:

```bibtex
@software{TorcelloRequena2025,
  author = {Torcello-Requena, Alberto},
  title = {MGall_Analysis: Comprehensive Bioinformatics Pipeline for 
           Mycoplasma gallisepticum Genomic Analysis},
  year = {2025},
  url = {https://github.com/[your-username]/MGall_Analysis}
}
```

**Key publications**:

- Phylogenetic methods: [IQ-TREE](http://www.iqtree.org/), [BEAST2](https://www.beast2.org/)
- Pangenome analysis: [Panaroo](https://github.com/gtonkinhill/panaroo)
- Recombination detection: [Gubbins](https://github.com/nickjcroucher/gubbins)

---

## License

This code is provided "as is" for research purposes. No warranty implied.

---

## Version History

- **v2.0** (2025-01-11): Complete refactoring with R integration
- **v1.0** (2024): Initial Python-only pipeline

---

**Last updated**: 2025-01-11
