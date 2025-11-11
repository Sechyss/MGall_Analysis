# MGall_Analysis

Comprehensive bioinformatics pipeline for *Mycoplasma gallisepticum* genomic analysis, including genome assembly, annotation, phylogenetics, pangenome analysis, GWAS, and epidemiological modeling.

**Author**: Dr. Alberto Torcello-Requena  
**Affiliation**: Environment and Sustainability Institute, University of Exeter  
**Contact**: alberto.torcello-requena@exeter.ac.uk

---

## Table of Contents

1. Overview
2. Repository Structure
3. Installation & Requirements
4. Workflow Overview
5. Detailed Script Documentation
6. R Analysis Pipeline
7. Typical Use Cases
8. Troubleshooting
9. Citation

---

## Overview

This repository contains a complete computational pipeline for analyzing *Mycoplasma gallisepticum* genomic data from house finches (*Haemorhous mexicanus*). The analysis spans:

- **Data management**: Sample manifest generation, metadata extraction, MD5 validation
- **Genome assembly**: De novo assembly (SPAdes), quality control
- **Annotation**: Prokka-based structural annotation
- **Phylogenetics**: Maximum likelihood trees (IQ-TREE), time-calibrated phylogenies (BEAST), recombination detection (Gubbins)
- **Pangenomics**: Gene presence/absence analysis (Panaroo), lineage-specific genes
- **Population genetics**: SNP calling, GWAS, Manhattan plots
- **Virulence analysis**: Phage detection (VIBRANT), defense systems (PADLOC)
- **Epidemiological modeling**: SEIRS models with parameter sweeps and sensitivity analysis

---

## Repository Structure

```
MGall_Analysis/
├── .idea/                          # PyCharm project settings
├── .vscode/                        # VS Code settings
├── MGall_R_analysis/               # R-based analyses
│   ├── 2_Functions/                # R utility functions
│   │   ├── 2_2_Lineage_detection_20240909.R
│   │   ├── 2_3_Lineage_fitness_20240909.R
│   │   └── 2_4_Lineage_defining_mutations.R
│   ├── Brownian_Alberto/           # Brownian motion trait evolution
│   │   ├── Brownian_Alberto.R
│   │   └── meta.txt                # Sample metadata
│   ├── BacDating.R                 # Bayesian dating analysis
│   ├── Panstripe.R                 # Pangenome stripe plots
│   ├── Phylowave.R                 # Phylogenetic wave analysis
│   ├── Position_segregating_sites.R
│   ├── Skyline_lineages.R          # Skyline plots per lineage
│   └── Skylineplot.R               # BEAST skyline plotting
├── Scripts/                        # Python & shell scripts
│   ├── General_python_scripts/     # General utilities
│   │   ├── Alignment_trimming.py
│   │   ├── Alignment_trimming_informativesites.py
│   │   ├── Alignment2Fasta.py
│   │   ├── AT_ratio_overtime.py
│   │   ├── GeneContext_verification.py
│   │   ├── GFFs2FNA.py
│   │   ├── Pangenome_translation.py
│   │   ├── Position_segregating_sites.py
│   │   ├── Preparation_trees.py
│   │   ├── Rename_dictionary.py
│   │   ├── Rename_files.py
│   │   ├── Rename_GFFs_metadata.py
│   │   ├── Rename_taxa_columns.py
│   │   ├── Rename_taxa_folders.py
│   │   └── Tree_VCF.py
│   ├── MG_genome/                  # Core genomic analyses
│   │   ├── Alignment_lineages.py
│   │   ├── Ancestral_tree_metadata.py
│   │   ├── Creation_list_input_gubbins.py
│   │   ├── Dictionary_clusters_heatmap.py
│   │   ├── GWAS_results_fasta_results.py
│   │   ├── Lineage_differences.py
│   │   ├── Lineage_differences_vcf.py
│   │   ├── Manhattan_annotated.py
│   │   ├── Manhattan_plot.py
│   │   ├── Manhattan_prep.py
│   │   ├── Manifest_files.py
│   │   ├── Metadata_extraction.py
│   │   ├── Metadata_Viral_analysis.py
│   │   ├── Pangenome_clustering_heatmap.py
│   │   ├── Pastml_analysis.py
│   │   ├── Phylogeny_snps.py
│   │   ├── Presence_absence_matrix.py
│   │   ├── SkyGrid_plot.py
│   │   └── SkyGrid_tree.py
│   ├── MG_motility/                # Motility gene analysis
│   │   ├── Analysis_vcf_motility.py
│   │   └── Motility_protein_preparation.py
│   ├── Shell_scripts/              # Bash pipeline scripts
│   │   ├── BCFtools.sh             # Consensus calling
│   │   ├── Download_sra.sh         # SRA data retrieval
│   │   ├── Gubbins_pipeline.sh     # Recombination detection
│   │   ├── Lucy_Prokka.sh          # Annotation pipeline
│   │   ├── Lucy_SPADES_Assembly.sh # De novo assembly
│   │   ├── Padloc.sh               # Defense system detection
│   │   ├── Panaroo_pipeline.sh     # Pangenome analysis
│   │   ├── Poppunk_code.sh         # PopPUNK clustering
│   │   ├── Prokka_command.sh
│   │   ├── Prokka_complete_genomes.sh
│   │   ├── SH_comparison.sh        # Tree topology tests
│   │   ├── Sickle_readTrimming.sh  # Quality trimming
│   │   ├── SPADES_Pipeline.sh
│   │   ├── Tree_formation.sh
│   │   ├── Treemaking_commands.sh
│   │   └── Vibrant.sh              # Phage prediction
│   └── Slurm_scripts/              # HPC cluster scripts
│       ├── IQ-Tree_slurm.sh
│       ├── Mgall_mapping_slurm.sh
│       └── Variant_calling.sh
└── README.md                       # This file
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

#### **Rename_taxa_columns.py** / **Rename_taxa_folders.py**
Batch renaming utilities for dataframes and file systems.

```bash
python Rename_taxa_columns.py --file input.tsv
python Rename_taxa_folders.py --directory /path/to/data/
```

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

---

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
```
Error: MD5 file contains extra text
```
**Solution**: Ensure `.md5` files contain only the 32-character hash:
```bash
md5sum file.fastq > file.fastq.md5
```

#### **2. GFF Annotation Failures**
```
Unknown gene (0bp away)
```
**Solutions**:
- Check chromosome name matching: `grep "^>" reference.fna` vs `awk '{print $1}' reference.gff`
- Verify coordinate system (GFF3 is 1-based)
- Increase distance window: `--distance 10000`

#### **3. Biopython Deprecation Warnings**
```
BiopythonDeprecationWarning: feature.strand is deprecated
```
**Solution**: Use `feature.location.strand` instead.

#### **4. Qt Backend Errors**
```
qt.qpa.plugin: Could not load the Qt platform plugin "xcb"
```
**Solution**: Set matplotlib backend in script header:
```python
import matplotlib
matplotlib.use('Agg')
```

#### **5. Prokka Fails on Large Genomes**
```
Error: Too many contigs
```
**Solution**: Filter short contigs before annotation:
```bash
seqkit seq -m 500 assembly.fasta > filtered.fasta
```

#### **6. IQ-TREE Memory Issues**
```
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