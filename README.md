# snp-calling-scripts

## Steps

## 02 Read Trimming
sbatch 02_read_trimming.sh setup/02_read_trimming.sh

## 03 Mapping using BWA-mem
03_run_mapping.py -r [PATH/TO/REFERENCE] -i [PATH/TO/TRIMMED/FASTQ/FILES] -o [PATH/TO/OUTPUT/FOLDERS] -p [PATH/TO/THE/SCRIPTS/FOLDER]