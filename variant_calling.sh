#!/bin/bash

#SBATCH --job-name="BQSR"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
#SBATCH --mem=25G
#SBATCH --partition=all
#SBATCH --output=BQSR-%j.out
#SBATCH --error=BQSR-%j.error


#create two vcf files of known_variants (SNPs and INDELs) with GATK [HaplotypeCaller] from _RG_sorted_dedup.bam file
#these two vcf files will be required for the base quality score recalibration process


#define attributes. Do not introduce extensions (e.g. .fastq.gz or .fa) in the sbatch command
REF=$1
SAMPLE=$2
BASEDIR=$3

#load required module
export PATH=/software/bin:$PATH;
module use /software/module/
module add vital-it
module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/samtools/1.8
module add UHTS/Analysis/picard-tools/2.18.11
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module add R/latest

#create a general working directory for 2_base_recalibration and go there
#create a working directory for individual {SAMPLE}_recalibration.bam files and go there
#create a local temporary directory for GATK4
mkdir -p ${BASEDIR}/2_base_recalibration
cd ${BASEDIR}/2_base_recalibration
mkdir base_recalibration_${SAMPLE}
cd base_recalibration_${SAMPLE}
mkdir temp

#link REF and _RG_sorted_dedup.bam files locally
ln -s ${REF} ${REF}.fa
ln -s /data/users/${USER}/1_RG_BAM_dedup/RG_BAM_dedup_${READS}/${READS}_RG_sorted_dedup.bam ${READS}_RG_sorted_dedup.bam

#Generate all the reference files required by GATK to work properly
REF="$(echo ${REF} | perl -pe 's/(\.\w+$)//')"
bwa index ${REF}.fa
samtools faidx ${REF}.fa
picard-tools CreateSequenceDictionary R=${REF}.fa O=${REF}.dict
picard-tools BuildBamIndex INPUT=${READS}_RG_sorted_dedup.bam

#create a shortcut for calling the GATK4 package and place it in a temporary directory
GATK4="java -Xmx25g -Djava.oi.tmpdir=`pwd`/tmp -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.0.4.0/bin/GenomeAnalysisTK.jar "

#generate a first list of raw SNPs and INDELs (${READS}_raw_variants_1.vcf) from the _RG_sorted_dedup.bam file with GATK [HaplotypeCaller]
$GATK4 HaplotypeCaller -R ${REF}.fa -I ${READS}_RG_sorted_dedup.bam -O ${READS}_raw_variants_1.vcf

#separate raw SNPs and INDELs in two separated vcf files since filtration parameters will be different for each (see further)
$GATK4 SelectVariants -R ${REF}.fa -V ${READS}_raw_variants_1.vcf -select-type SNP -O ${READS}_raw_SNPs_1.vcf
$GATK4 SelectVariants -R ${REF}.fa -V ${READS}_raw_variants_1.vcf -select-type INDEL -O ${READS}_raw_INDELs_1.vcf

#hard-filter raw SNPs and INDELs vcf files with [VariantFiltration] to mark the bad ones with FILTER (many variants in both raw files are not genuine variants)
#SNPs and INDELs matching any of the specified criteria will be considered bad and marked FILTER
#for _raw_SNPs_1
$GATK4 VariantFiltration -R ${REF}.fa -V ${READS}_raw_SNPs_1.vcf -O ${READS}_raw_SNPs_2.vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 35.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name "LMS_SNP_filter1"
#for _raw_ INDELs_1
$GATK4 VariantFiltration -R ${REF}.fa -V ${READS}_raw_INDELs_1.vcf -O ${READS}_raw_INDELs_2.vcf --filter-expression 'QD < 2.0 || FS > 200.0' --filter-name "LMS_INDEL_filter1"


#Part 4: recalibration of base quality score with filtered {READS}_raw_SNPs_2.vcf and {READS}_raw_INDELs_2.vcf


#iteration 1
#use {READS}_raw_SNPs_2.vcf and {READS}_raw_INDELs_2.vcf to generate a first model of recalibration (tables;_recal_model_1.grp) that will be used for base quality score recalibration
$GATK4 BaseRecalibrator -R ${REF}.fa -I ${READS}_RG_sorted_dedup.bam --known-sites ${READS}_raw_SNPs_2.vcf --known-sites ${READS}_raw_INDELs_2.vcf -O ${READS}_recal_model_1.grp

#recalibrate the original bam file (GATK4 does not allow the recalibration on-the-fly like GATK3)
$GATK4 ApplyBQSR -R ${REF}.fa -I ${READS}_RG_sorted_dedup.bam --bqsr-recal-file ${READS}_recal_model_1.grp -O ${READS}_recal_1.bam

#use the recal_1.bam to built de second model
$GATK4 BaseRecalibrator -R ${REF}.fa -I ${READS}_recal_1.bam --known-sites ${READS}_raw_SNPs_2.vcf --known-sites ${READS}_raw_INDELs_2.vcf -O ${READS}_recal_model_2.grp

#plots the two models of recalibration to see whether convergence is observed or not
$GATK4 AnalyzeCovariates -before ${READS}_recal_model_1.grp -after ${READS}_recal_model_2.grp -plots ${READS}_recalibration_1_2.pdf


#Part 5 : Call variants on the recalibrated bam file

#this generates a ${READS}_final_variants.g.vcf file that contains both SNPs and INDELs
#runs HC with the -ERC GVCF mode to add relevant informations for downstream genotypingt
$GATK4 HaplotypeCaller -R ${REF}.fa -I ${READS}_recal_1.bam -O ${READS}_final_variants.g.vcf -ERC GVCF --pcr-indel-model NONE

