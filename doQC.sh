#!/bin/bash

#SBATCH --job-name=QC # Job name
#SBATCH --time=2-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=all # Partition to submit to
#SBATCH --output=/home/ubelix/ips/rchoudhury/Data/biscut_ind/trimlog/QC-\%A_\%a.out # File to which STDOUT will be written
#SBATCH --error=/home/ubelix/ips/rchoudhury/Data/biscut_ind/trimlog/QC-\%A_\%a.err # File to which STDERR will be written 
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=rimjhim.choudhury@ips.unibe.ch # Email to which notifications will be sent
#SBATCH --array=1-31                 # Array range

# - Rimjhim Roy Choudhury 
# This is an script that combines array tasks with
# bash loops to process many short runs. Array jobs are convenient
# for running lots of tasks, but if each task is short, they
# quickly become inefficient, taking more time to schedule than
# they spend doing any work and bogging down the scheduler for
# all users.
 
pwd; hostname; date
start=`date +%s`

#Set the number of runs that each SLURM task should do
PER_TASK=6


echo "Hostname: $(eval hostname)"


if [ -z ${SLURM_NTASKS} ]; then THREADS=$( nproc ); else THREADS=${SLURM_NTASKS}; fi

if [ -z $1 ]; then echo "json path is unset"; exit 1; else echo "json setup file: '$1'"; fi

module add vital-it UHTS/Analysis/trimmomatic/0.36
#source /data/users/lfalquet/BC7107_18/scripts/module.sh
#source ~/.bash_profile
source activate py3
#input
basepath=$(eval jq .project.basename $1 | sed 's/^"\(.*\)"$/\1/')
id=$(eval jq .project.id $1 | sed 's/^"\(.*\)"$/\1/' )
fastqs=$(eval jq .project.fastqs $1 | sed 's/^"\(.*\)"$/\1/')
samples=$(eval jq .project.samples $1 | sed 's/^"\(.*\)"$/\1/')
#output
qc_out="${basepath}qc_trim"

mkdir $qc_out
cd $qc_out
echo sample is $samples

declare -a sam=(`cat $samples`)
#sam=(`cat $sample`)


# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))
 
# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM
 
# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  sam=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $samples)

  trun="trimmomatic PE -phred33 -threads 12 ${fastqs}/${sam[$run]}_R1_001.fastq.gz ${fastqs}/${sam[$run]}_R2_001.fastq.gz ${qc_out}/${sam[$run]}_R1_001.trim.fastq.gz ${qc_out}/${sam[$run]}_R1_001.unpaired.fastq.gz ${qc_out}/${sam[$run]}_R2_001.trim.fastq.gz ${qc_out}/${sam[$run]}_R2_001.unpaired.fastq.gz ILLUMINACLIP:$HOME/Data/adapters/TruSeq.all.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"
echo $trun
eval $trun
done
 

## Commented script follows:
: <<'END_COMMENT'
for i in `ls -1 ${fastqs}/*.fastq.gz`
do
call1="fastqc -t 4 ${i}"
echo $call1
eval $call1
done
END_COMMENT

: <<'END_COMMENT'
for k in $(echo $samples| sed "s/,/ /g")
do
call2="/home/rchoudhury/programs/FastQC/fastqc -t 4 ${qc_out}/${k}_*.trim.fastq.gz"
echo $call2
eval $call2
done

multi="multiqc ."
echo $multi
eval $multi
END_COMMENT






