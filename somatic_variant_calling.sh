#!/bin/bash

#set working directory
# cd /mnt/d/courses/genomics_tut/somatic_vatriant_call

# Create directories if required
# mkdir data
# mkdir reference
# mkdir aligned
# mkdir results
# mkdir mutect2

# Sebset the reads (This is for practice purpose on you local machines as entire datasets cannot be run on local machines with basic ram)

# seqtk sample -s100 SLGFSK-N_231335_r1_chr5_12_17.fastq.gz 100000 > data/normal_R1.fastq
# seqtk sample -s100 SLGFSK-N_231335_r2_chr5_12_17.fastq.gz 100000 > data/normal_R2.fastq
# seqtk sample -s100 SLGFSK-T_231336_r1_chr5_12_17.fastq.gz 100000 > data/tumor_R1.fastq
# seqtk sample -s100 SLGFSK-T_231336_r2_chr5_12_17.fastq.gz 100000 > data/tumor_R2.fastq


# GATK4 best practices

#--------------------------------------------------------------------------#

# Files to be downloaded only once (human genome reference, indexing, dictionary generation, file for known sites for BQSR, mutect2 supporting files)

# Human genome reference file download and unzip
# wget -P reference/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# gunzip reference/hg38.fa.gz

# Index the reference file
# samtools faidx reference/hg38.fa

# reference dict.: required before running haplotype caller
# gatk CreateSequenceDictionary R= reference/hg38.fa O= reference/hg38.dict

# Download known sites files for BQSR from GATK resource bundle
# Download dbSNP VCF
#wget -P ~/Desktop/demo/reference/ https://storage.googleapis.com/gatk-best-practices/somatic-hg38/Homo_sapiens_assembly38.dbsnp138.vcf

# Download dbSNP index
#wget -P ~/Desktop/demo/reference/ https://storage.googleapis.com/gatk-best-practices/somatic-hg38/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Mutect2 supporting files

# gnomAD
# wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz ~/Desktop/demo/supporting_files/mutect2_supporting_files
# wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi ~/Desktop/demo/supporting_files/mutect2_supporting_files

# PoN
# wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz ~/Desktop/demo/supporting_files/mutect2_supporting_files
# wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi ~/Desktop/demo/supporting_files/mutect2_supporting_files

# For creating panel of normals
# wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi ~/Desktop/demo/supporting_files/mutect2_supporting_files

# Intervals list
# wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list ~/Desktop/demo/supporting_files/mutect2_supporting_files



#--------------------------------------------------------------------------#

# Set directories path tp variable so that the it will be easy to include them in the commands
ref=/mnt/d/courses/genomics_tut/somatic_vatriant_call/reference/hg38.fa
reference=/mnt/d/courses/genomics_tut/somatic_vatriant_call/reference
known_sites=/mnt/d/courses/genomics_tut/somatic_vatriant_call/reference/Homo_sapiens_assembly38.dbsnp138.vcf
proj_dir=/mnt/d/courses/genomics_tut/somatic_vatriant_call
aligned=$proj_dir/aligned
reads=$proj_dir/data
results=$proj_dir/results
gatk_path=/mnt/d/courses/genomics_tut/gatk-4.5.0.0/gatk
mutect2_files=/mnt/d/courses/genomics_tut/somatic_vatriant_call/mutect2_files
data_source_path=/mnt/d/courses/genomics_tut/ensembl_vep/

#--------------------------------------------------------------------------#
# Preprocessing files

# ---1. QC ---
 fastqc ${reads}/normal_R1.fastq -o ${reads}
 fastqc ${reads}/normal_R2.fastq -o ${reads}
 fastqc ${reads}/tumor_R1.fastq -o ${reads}
 fastqc ${reads}/tumor_R2.fastq -o ${reads}

# ---2. Alignment ---

# BWA index reference
 bwa index ${ref}

# BWA alignment
 bwa mem -t 4 -R "@RG\tID:HT33CBBXX.3\tPL:ILLUMINA\tSM:normal" ${ref} ${reads}/normal_R1.fq ${reads}/normal_R2.fast > ${aligned}/normal.paired.sam
 bwa mem -t 4 -R "@RG\tID:HT33CBBXX.3\tPL:ILLUMINA\tSM:tumor" ${ref} ${reads}/tumor_R1.fq ${reads}/tumor_R2.fast > ${aligned}/tumor.paired.sam

 bwa mem -t 4 -R "@RG\tID:HT33CBBXX.3\tPL:ILLUMINA\tSM:normal" ${ref} ${reads}/normal_R1.fastq ${reads}/normal_R2.fastq > ${aligned}/normal.paired.sam
 bwa mem -t 4 -R "@RG\tID:HT33CBBXX.3\tPL:ILLUMINA\tSM:tumor" ${ref} ${reads}/tumor_R1.fastq ${reads}/tumor_R2.fastq > ${aligned}/tumor.paired.sam

# --- Mark duplicates and sort ---
 ${gatk} MarkDuplicatesSpark -I ${aligned}/normal.paired.sam -O ${aligned}/normal_sorted_dedup.bam
 ${gatk} MarkDuplicatesSpark -I ${aligned}/tumor.paired.sam -O ${aligned}/tumor_sorted_dedup.bam

# --- Base quality recaliberation ---
# Build model
 ${gatk} BaseRecalibrator -I ${aligned}/normal_sorted_dedup.bam -R ${ref} --known-sites ${reference}/Homo_sapiens_assembly38.dbsnp138.vcf -O ${aligned}/normal_recal.table
 ${gatk} BaseRecalibrator -I ${aligned}/tumor_sorted_dedup.bam -R ${ref} --known-sites ${reference}/Homo_sapiens_assembly38.dbsnp138.vcf -O ${aligned}/tumor_recal.table

# --- Apply model to adjust base quality scores ---
 ${gatk} ApplyBQSR -I ${aligned}/normal_sorted_dedup.bam -R ${ref} --bqsr-recal-file ${aligned}/normal_recal.table -O ${aligned}/normal_sorted_dedup_bqsr.bam
 ${gatk} ApplyBQSR -I ${aligned}/tumor_sorted_dedup.bam -R ${ref} --bqsr-recal-file ${aligned}/tumor_recal.table -O ${aligned}/tumor_sorted_dedup_bqsr.bam

# -- Collect alignment and insert size metric ---
# These steps are done for quality control and validation of your sequencing/alignment before variant calling
# Install R for histogram
 ${gatk} CollectAlignmentSummaryMetrics R=${ref} I=${aligned}/normal_sorted_dedup_bqsr.bam O=${aligned}/normal_alignment_metrics.txt
 ${gatk} CollectInsertSizeMetrics INPUT=${aligned}/normal_sorted_dedup_bqsr.bam OUTPUT=${aligned}/normal_insert_size_metrics.txt HISTOGRAM_FILE=${aligned}/normal_insert_size_histogram.pdf

 ${gatk} CollectAlignmentSummaryMetrics R=${ref} I=${aligned}/tumor_sorted_dedup_bqsr.bam O=${aligned}/tumor_alignment_metrics.txt
 ${gatk} CollectInsertSizeMetrics INPUT=${aligned}/tumor_sorted_dedup_bqsr.bam OUTPUT=${aligned}/tumor_insert_size_metrics.txt HISTOGRAM_FILE=${aligned}/tumor_insert_size_histogram.pdf

#--------------------------------------------------------------------------#

# Somatic Variant Calling

# Call variants
 ${gatk} gatk Mutect2 -R ${ref} -I ${aligned}/tumor_sorted_dedup_bqsr.bam -I ${aligned}/normal_sorted_dedup_bqsr.bam -tumor tumor -normal normal --germline-resource ${mutect2_files}/af-only-gnomad.hg38.vcf.gz --panel-of-normals ${mutect2_files}/1000g_pon.hg38.vcf.gz -O $results/somatic_vatiants_mutect2.vcf.gz --f1r2-tar-gz ${results}/variants_f1r2.tar.gz

# Cross sample contamination estimation
# GetPileupSummaries: Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results are used with CalculateContamination.

# For tumor
 ${gatk_path} GetPileupSummaries --java-options '-Xmx50G' -I ${aligned}/tumor_sorted_dedup_bqsr.bam -V ${mutect2_files}/af-only-gnomad.hg38.vcf.gz -L ${mutect2_files}/exome_calling_regions.v1.1.interval_list -O ${results}/tumor_getpilesummaries.table
# For normal
 ${gatk_path} GetPileupSummaries --java-options '-Xmx50G' -I ${aligned}/normal_sorted_dedup_bqsr.bam -V ${mutect2_files}/af-only-gnomad.hg38.vcf.gz -L ${mutect2_files}/exome_calling_regions.v1.1.interval_list -O ${results}/normal_getpilesummaries.table

# Calculate contamination
 ${gatk_path} CalculateContamination -I ${results}/tumor_getpilesummaries.table -matched ${results}/normal_getpilesummaries.table -O ${results}/pair_calculatecontamination.table

# Estimate read orientation artifacts
 ${gatk_path} LearnReadOrientationModel -I ${results}/variants_f1r2.tar.gz -O ${results}/read-orientation-model.tar.gz

# Filtering variants
 ${gatk_path} FilterMutectCalls -V ${results}/somatic_vatiants_mutect2.vcf.gz -R ${ref} --contamination-table ${results}/pair_calculatecontamination.table --ob-priors ${results}/read-orientation-model.tar.gz -O ${results}/filtered_somatic_variants_mutect2.vcf

# Variant annotation by funcotator
# Download funcotator data sources from google cloud which requires billing account
# ${gatk_path} Funcotator --variant ${results}/filtered_somatic_variants_mutect2.vcf --reference ${ref} --ref-version hg38 --data-source-path ${data_source_path} --output ${results}/somatic_variants_funcotated.vcf --output-file-format VCF

# ## OR ELSE ##

# Variant annotation by Ensemble VEP
# funcotator data sources is now availble from google cloud which requires a billing account, hence we are using Ensembl VEP
# Download ensembl vep  data source, follow the steps below 
# conda install -c bioconda ensembl-vep
# vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/vep_data

# vep -i ${results}/filtered_somatic_variants_mutect2.vcf -o ${results} annotated.vep.vcf --cache --dir_cache ${data_source_path} --assembly GRCh38 --vcf