# Rare Genetic Disease Diagnosis in Family Trios – Genomic Analysis Pipeline

## Authors

Laura Mottarlini & Sabrina Pervin Abu

## Introduction

This project presents a complete bioinformatics pipeline for the diagnosis of rare Mendelian autosomal diseases within family trios (mother, father, child). The `script_genomics.sh` automates quality control, sequence alignment, variant calling, and filtering based on inheritance models, leveraging a combination of established tools in a streamlined workflow.

Our aim is to distinguish between de novo and inherited variants and evaluate their pathogenic potential based on genomic annotation and impact scoring.

## Overview

The `script_genomics.sh` pipeline includes the following main stages:

1. **Input Preparation**:
   - Defines a set of case IDs and sample roles (mother, father, child).
   - Verifies and organizes input FASTQ files.

2. **Quality Control**:
   - Uses `FastQC` and `Qualimap` to assess sequence and alignment quality.
   - Aggregates reports with `MultiQC`.

3. **Sequence Alignment**:
   - Aligns reads to the reference genome using `Bowtie2`.
   - Converts and sorts alignment outputs via `samtools`.

4. **Coverage Analysis**:
   - Computes BedGraph coverage with `bedtools genomecov`.

5. **Variant Calling**:
   - Generates trio-based multi-sample VCF files using `FreeBayes`.
   - Applies thresholds for mapping and base quality.

6. **Variant Prioritization**:
   - Filters variants using specific genotype patterns:
     - AD: `0/1` (child), `0/0` (parents)
     - AR: `1/1` (child), `0/1` (parents)
   - Intersects variants with exonic regions using `bedtools`.

7. **Annotation & Visualization**:
   - Annotates VCFs with VEP (Ensembl, GRCh37).
   - Visualizes variants on the UCSC Genome Browser.

## Usage

Run the script with:

```bash
chmod +x script_genomics.sh
nohup ./script_genomics.sh &
```

## Script Workflow

The pipeline performs the following steps for each case:

1. **FASTQ Quality Check**  
   - Tool: `FastQC`  
   - Assesses the quality of raw sequencing reads.

2. **Alignment to Reference Genome**  
   - Tool: `Bowtie2`  
   - Aligns reads from mother, father, and child to the reference genome.

3. **SAM to Sorted BAM Conversion**  
   - Tool: `samtools`  
   - Converts alignment output from SAM to BAM, and sorts the files.

4. **Coverage and Exon Intersection**  
   - Tool: `bedtools genomecov`, `bedtools intersect`  
   - Computes read coverage and intersects variants with exonic regions.

5. **Variant Calling**  
   - Tool: `FreeBayes`  
   - Generates multi-sample VCFs based on trio alignments.

6. **Variant Filtering**  
   - Tool: `grep`  
   - Filters variants based on specific genotype patterns:
     - Autosomal Dominant (AD): `0/1` in child, `0/0` in parents.
     - Autosomal Recessive (AR): `1/1` in child, `0/1` in parents.

7. **Multi-Sample VCF Sorting**  
   - Tool: `bcftools`  
   - Ensures consistent sample order in VCF headers.

8. **Report Aggregation**  
   - Tool: `MultiQC`  
   - Merges quality control results into a unified report.

9. **Manual Annotation Review**  
   - Tools: `VEP` (Ensembl Variant Effect Predictor), `UCSC Genome Browser`  
   - Used to assess variant impact, frequency, and clinical relevance.

---

## Dependencies

Make sure the following tools are installed and accessible in your environment:

- `bowtie2`
- `samtools`
- `freebayes`
- `bcftools`
- `bedtools`
- `fastqc`
- `qualimap`
- `multiqc`
- `vep` *(for annotation, run externally)*

---

## Input File Format

The script expects input files for each case to be named as follows:

- `case<caseID>_mother.fq.gz`
- `case<caseID>_father.fq.gz`
- `case<caseID>_child.fq.gz`
