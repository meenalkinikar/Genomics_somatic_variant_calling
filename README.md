# Somatic Variant Calling Pipeline

A reproducible somatic variant calling workflow for paired tumor-normal sequencing data following GATK Best Practices.

## Tools Used

- FastQC
- Seqtk
- BWA-MEM
- SAMtools
- Picard
- GATK 4 (Mutect2, BQSR, Funcotator)

## Workflow

1. Read quality assessment using FastQC
2. Read subsetting for testing using Seqtk
3. Alignment to the human reference genome (hg38) using BWA-MEM
4. BAM sorting and indexing
5. Duplicate marking using Picard
6. Base Quality Score Recalibration (BQSR)
7. Somatic variant calling using Mutect2
8. Variant filtering and contamination assessment
9. Functional annotation using Funcotator/Ensemble VEP

## Input

- Paired-end FASTQ files from matched normal and tumor samples
- Human reference genome (hg38)
- Known variant resources (dbSNP, gnomAD)

## Output

- Recalibrated BAM files
- Somatic SNV/Indel VCF files
- Filtered somatic variants
- Functionally annotated variants


