#!/bin/bash

#SBATCH --job-name=QC # Job name
#SBATCH --nodes=1
#SBATCH --time=2-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=all # Partition to submit to
#SBATCH --output=QC-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=QC-%N-%j.err # File to which STDERR will be written #SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=rimjhim.choudhury@ips.unibe.ch # Email to which notifications will be sent

#create and go to the TP directory
start=`date +%s`



echo "Hostname: $(eval hostname)"

if [ -z ${SLURM_JOB_ID} ]; then echo "running locally"; else echo "My SLURM_JOB_ID: ${SLURM_JOB_ID}"; fi
if [ -z ${SLURM_NTASKS} ]; then THREADS=$( nproc ); else THREADS=${SLURM_NTASKS}; fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi


module add vital-it UHTS/Quality_control/fastqc
conda activate defaults3
#input
basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
fastqs=$(eval jq .project.fastqs $1 | sed 's/^"\(.*\)"$/\1/')

#output
qc_out="${basepath}/fastqc"

mkdir $qc_out
cd $qc_out



for i in `ls -1 ${fastqs}/*.fastq.gz`
do
call1="fastqc -t 4 ${i}"
echo $call1
eval $call1
done


multi="multiqc ."
echo $multi
eval $multi