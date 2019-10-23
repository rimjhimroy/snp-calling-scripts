#!/bin/bash

#SBATCH --job-name="join-genotyping"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=92:00:00
#SBATCH --mem=48G
#SBATCH --partition=all
#SBATCH --error=join_genotyping-%j.error


#consolidate several g.vcf files before the join genotyping step
#the sbatch command will be sbatch scriptC-bioinfo.slurm ATH READS1 READS2 READS3 ... READS8 GENOTYPING_GROUP
#exemple of GENOTYPING_GROUP: vitE
REF=$1
READS1=$2
READS2=$3
READS3=$4
READS4=$5
READS5=$6
READS6=$7
READS7=$8
READS8=$9
GENOTYPING_GROUP=${10}
THREADS=8


#load required modules
export PATH=/software/bin:$PATH;
module use /software/module/
module add vital-it
module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/samtools/1.4
module add UHTS/Analysis/picard-tools/2.9.0
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0


#create a general working directory for genotyping  and go there. -p: create the directory if it does not exist
#create a second directory with the genotyping project name and go there
#create a local temporary directory for GATK4
mkdir -p /data/users/${USER}/3_genotyping
cd /data/users/${USER}/3_genotyping
mkdir -p ${GENOTYPING_GROUP}
cd ${GENOTYPING_GROUP}
mkdir -p temp


#link REF and ${READS}_final_variants.g.vcf files locally
#don't forget the index files that are required
ln -s /data/users/lfalquet/BC7107_18/reference/${REF}.fa ${REF}.fa
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS1}/${READS1}_final_variants.g.vcf ${READS1}_final_variants.g.vcf
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS2}/${READS2}_final_variants.g.vcf ${READS2}_final_variants.g.vcf
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS3}/${READS3}_final_variants.g.vcf ${READS3}_final_variants.g.vcf
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS4}/${READS4}_final_variants.g.vcf ${READS4}_final_variants.g.vcf
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS5}/${READS5}_final_variants.g.vcf ${READS5}_final_variants.g.vcf
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS6}/${READS6}_final_variants.g.vcf ${READS6}_final_variants.g.vcf
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS7}/${READS7}_final_variants.g.vcf ${READS7}_final_variants.g.vcf
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS8}/${READS8}_final_variants.g.vcf ${READS8}_final_variants.g.vcf

ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS1}/${READS1}_final_variants.g.vcf.idx ${READS1}_final_variants.g.vcf.idx
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS2}/${READS2}_final_variants.g.vcf.idx ${READS2}_final_variants.g.vcf.idx
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS3}/${READS3}_final_variants.g.vcf.idx ${READS3}_final_variants.g.vcf.idx
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS4}/${READS4}_final_variants.g.vcf.idx ${READS4}_final_variants.g.vcf.idx
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS5}/${READS5}_final_variants.g.vcf.idx ${READS5}_final_variants.g.vcf.idx
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS6}/${READS6}_final_variants.g.vcf.idx ${READS6}_final_variants.g.vcf.idx
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS7}/${READS7}_final_variants.g.vcf.idx ${READS7}_final_variants.g.vcf.idx
ln -s /data/users/${USER}/2_base_recalibration/base_recalibration_${READS8}/${READS8}_final_variants.g.vcf.idx ${READS8}_final_variants.g.vcf.idx


REF="$(echo ${REF} | perl -pe 's/(\.\w+$)//')"
bwa index ${REF}.fa
samtools faidx ${REF}.fa
picard-tools CreateSequenceDictionary R=${REF}.fa O=${REF}.dict

#create a shortcut for calling the GATK4 package and place it in a temporary directory
GATK4="java -Xmx46g -Xms46g -Djava.oi.tmpdir=`pwd`/tmp -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.0.4.0/bin/GenomeAnalysisTK.jar "


#Merge the 8 g.vcf files to genotype them with [GenomicsBDImport]
#at the moment, this can be done only for one interval [chromosome]. "For now you need to run on each interval separately"

#chromosome 1
$GATK4 GenomicsDBImport -V ${READS1}_final_variants.g.vcf -V ${READS2}_final_variants.g.vcf -V ${READS3}_final_variants.g.vcf -V ${READS4}_final_variants.g.vcf -V ${READS5}_final_variants.g.vcf -V ${READS6}_final_variants.g.vcf \
-V ${READS7}_final_variants.g.vcf -V ${READS8}_final_variants.g.vcf --reader-threads 8 --genomicsdb-workspace-path ${GENOTYPING_GROUP}_chr1_db -L 1

#chromosome 2
$GATK4 GenomicsDBImport -V ${READS1}_final_variants.g.vcf -V ${READS2}_final_variants.g.vcf -V ${READS3}_final_variants.g.vcf -V ${READS4}_final_variants.g.vcf -V ${READS5}_final_variants.g.vcf -V ${READS6}_final_variants.g.vcf \
-V ${READS7}_final_variants.g.vcf -V ${READS8}_final_variants.g.vcf --reader-threads 8 --genomicsdb-workspace-path ${GENOTYPING_GROUP}_chr2_db -L 2

#chromosome 3
$GATK4 GenomicsDBImport -V ${READS1}_final_variants.g.vcf -V ${READS2}_final_variants.g.vcf -V ${READS3}_final_variants.g.vcf -V ${READS4}_final_variants.g.vcf -V ${READS5}_final_variants.g.vcf -V ${READS6}_final_variants.g.vcf \
-V ${READS7}_final_variants.g.vcf -V ${READS8}_final_variants.g.vcf --reader-threads 8 --genomicsdb-workspace-path ${GENOTYPING_GROUP}_chr3_db -L 3

#chromosome 4
$GATK4 GenomicsDBImport -V ${READS1}_final_variants.g.vcf -V ${READS2}_final_variants.g.vcf -V ${READS3}_final_variants.g.vcf -V ${READS4}_final_variants.g.vcf -V ${READS5}_final_variants.g.vcf -V ${READS6}_final_variants.g.vcf \
-V ${READS7}_final_variants.g.vcf -V ${READS8}_final_variants.g.vcf --reader-threads 8 --genomicsdb-workspace-path ${GENOTYPING_GROUP}_chr4_db -L 4

#chromosome 5
$GATK4 GenomicsDBImport -V ${READS1}_final_variants.g.vcf -V ${READS2}_final_variants.g.vcf -V ${READS3}_final_variants.g.vcf -V ${READS4}_final_variants.g.vcf -V ${READS5}_final_variants.g.vcf -V ${READS6}_final_variants.g.vcf \
-V ${READS7}_final_variants.g.vcf -V ${READS8}_final_variants.g.vcf --reader-threads 8 --genomicsdb-workspace-path ${GENOTYPING_GROUP}_chr5_db -L 5


#Genotype merged gVCF files [GenotypeGVCFs] and group them all together [MergeVcfs]
#since these (files) databases were generated separately (per chromosome), I guess I have to run the command for each interval and merge them later (I guess)
$GATK4 GenotypeGVCFs -R ${REF}.fa -V gendb://${GENOTYPING_GROUP}_chr1_db -O ${GENOTYPING_GROUP}_chr1.vcf.gz
$GATK4 GenotypeGVCFs -R ${REF}.fa -V gendb://${GENOTYPING_GROUP}_chr2_db -O ${GENOTYPING_GROUP}_chr2.vcf.gz
$GATK4 GenotypeGVCFs -R ${REF}.fa -V gendb://${GENOTYPING_GROUP}_chr3_db -O ${GENOTYPING_GROUP}_chr3.vcf.gz
$GATK4 GenotypeGVCFs -R ${REF}.fa -V gendb://${GENOTYPING_GROUP}_chr4_db -O ${GENOTYPING_GROUP}_chr4.vcf.gz
$GATK4 GenotypeGVCFs -R ${REF}.fa -V gendb://${GENOTYPING_GROUP}_chr5_db -O ${GENOTYPING_GROUP}_chr5.vcf.gz


#goup individual ${GENOTYPING_GROUP}_chr*.vcf.gz files together
picard-tools MergeVcfs I=${GENOTYPING_GROUP}_chr1.vcf.gz I=${GENOTYPING_GROUP}_chr2.vcf.gz I=${GENOTYPING_GROUP}_chr3.vcf.gz I=${GENOTYPING_GROUP}_chr4.vcf.gz I=${GENOTYPING_GROUP}_chr5.vcf.gz O=${GENOTYPING_GROUP}_allchr.vcf.gz
