# CUT&RUN Analysis Pipeline

This repository contains a Snakemake pipeline for CUT&RUN data processing and analysis, including quality control, mapping, peak calling, and data visualization.

---

## 📖 Overview

This pipeline automates the following steps for CUT&RUN data analysis:

1. **Data Preparation and Quality Control:**
   - Adapter trimming using `trim_galore`
   - Quality control using `FastQC`

2. **Alignment and BAM Processing:**
   - Read alignment using `bowtie2`
   - Sorting and indexing with `samtools`
   - Quality assessment with `samtools idxstats`, `flagstat`, and `stats`
   - Duplicate marking and removal using `picard`

3. **Peak Calling and Signal Processing:**
   - Peak calling using `MACS3` and `SEACR`
   - BigWig file generation using `bedtools` and `bedGraphToBigWig`

4. **Data Visualization:**
   - Quality metrics visualization using `deeptools`
   - Fingerprint plots and PCA analysis for quality assessment

---

## 📦 Requirements

### Tools:
- Python (3.x)
- Snakemake
- Bowtie2
- Samtools
- Picard
- MACS3
- SEACR
- Bedtools
- DeepTools
- TrimGalore

### Python Packages:
- `pandas`
- `pysam`

**Ensure all tools are properly installed and available in your system path.**

---

## 📂 Folder Structure

```plaintext
├── config/
│   └── config.yaml         # Configuration file for the pipeline
├── rules/                  # Snakemake rules for each step
│   ├── common.smk
│   ├── trim.smk
│   ├── align.smk
│   ├── peak.smk
│   ├── qc.smk
│   └── deeptools.smk
├── data/                   # Raw sequencing data files (not included)
├── results/                # Output results directory
├── qc/                     # Quality control outputs
└── Snakefile               # Main entry point for the pipeline
```

---


##  📜 Input File
The pipeline requires a **sample metadata file** in TSV format. Below is an example:

| Name       | Unit | Fastq1                                        | Fastq2                                        | Control |
|------------|------|-----------------------------------------------|-----------------------------------------------|---------|
| MGG152 | KDM_KO_igg   | /groups/lackgrp/projects/col-cutandrun-tbo/rawdata/30-1110336138/00_fastq/C10_R1_001.fastq.gz | /groups/lackgrp/projects/col-cutandrun-tbo/rawdata/30-1110336138/00_fastq/C10_R2_001.fastq.gz | -   |
| MGG152 | KDM_KO_H3K27me3 | /groups/lackgrp/projects/col-cutandrun-tbo/rawdata/30-1110336138/00_fastq/C11_R1_001.fastq.gz | /groups/lackgrp/projects/col-cutandrun-tbo/rawdata/30-1110336138/00_fastq/C11_R2_001.fastq.gz | MGG152.KDM_KO_igg   |


- **Name**: Group name.
- **Unit**: Unit or replicate identifier.
- **Fastq1**: Path to raw FASTQ file (Forward Pair).
- **Fastq2**: Path to raw FASTQ file (Reverse Pair).
- **Control**: "\<Name\>.\<Unit\>" of the sample.
---

---

## 🎛️ Configurations
The configuration file (`config/config.yaml`) specifies pipeline parameters. Below is an example:

```yaml
SAMPLES: config/samples.tsv

OUTPUT:
    REF: hg38
    RUN:
        QC: True
        PEAKS: True
        BWS: True
    SAECR_MODE: 
        - stringent
        - relaxed


# Reference Genome Settings
REFERENCES:
    hg38:
        FA: /groups/lackgrp/genomeAnnotations/hg38/hg38.fa
        BWA_IDX: /groups/lackgrp/genomeAnnotations/hg38/hg38.bwa.idx
        BOWTIE2_IDX: /groups/lackgrp/genomeAnnotations/hg38/hg38.bowtie2.idx/hg38
        CHROM_SIZES: /groups/lackgrp/genomeAnnotations/hg38/hg38.chrom.sizes
        BLACKLIST: /groups/lackgrp/genomeAnnotations/hg38/hg38-blacklist.v2.bed
```

---



## ▶️ Running the Pipeline

### 1. Configure the Pipeline
Edit the `config/config.yaml` file to specify:
- Reference genome details.
- Output directory structure.
- Flags to enable or disable specific steps (e.g., QC, peak calling, BigWig normalization).

Ensure that the `config/samples.tsv` file is properly formatted with the sample information.

---

### 2. Slurm Profile
The pipeline uses a Slurm cluster via the `profile/` directory. The `config.yaml` for the Slurm profile should include the following:

```yaml
cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --partition={resources.partition}
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=smk-{rule}-{wildcards}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
default-resources:
  - partition=normal,big-mem,long,express
  - mem_mb=700000
  - disk_mb=1024000
restart-times: 1
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 12
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True

```
- **Cluster Resources**: Adjust memory (mem_mb), disk (disk_mb), and partition names according to your Slurm setup.
- **Logging**: Logs for each rule are stored in logs/{rule}/.

---
### 3. Submission Script

Use the following run_pipeline.sh script to submit the pipeline to the Slurm cluster. The script activates the required conda environment and runs Snakemake with the specified profile.

```bash
#!/bin/bash
#SBATCH -c 64
#SBATCH --mem 720GB
#SBATCH -p long,big-mem,normal,express

source ~/.bashrc
conda activate cutandrun

snakemake --profile profile/

```
---
### 4. Submit the Pipeline

Run the following command to execute the pipeline:

```bash
sbatch run_pipeline.sh
```
This will:
- Automatically submit jobs to the Slurm cluster.
- Use the configuration specified in the profile/config.yaml file.
- Execute all defined rules in the pipeline.