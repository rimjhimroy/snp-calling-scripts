# Snp calling scripts

## Preprocess

You should not merge the different FASTQs from different lanes before aligning them. If you do so, then you will loose "read group" (RG tag) information. I would align each lane separately, then the scripts add specific RG tags to each lane and then merge the aligned bams (respecting RG information). RG information is important for downstream analyses in GATK, such as  BQSR.  

1. Rename [SAMPLENAME]_[LANE]_[PAIR]_[random_string].fastq.gz files by adding the run id to the file name, giving [SAMPLE]_[LANE]_[PAIR]_001.fastq.gz (important if there are several runs from the same library, else ignore). Make sure there is no '_' in SAMPLENAME.  

```
cd raw
for i in $( ls *.gz ); do x=$(zcat $i | head -n 1 | awk -F ":" '{print $2}');  sam=$(echo $i | awk -F "_" '{print $1}'); lane=$(echo $i | awk -F "_" '{print $2}');read=$(echo $i | awk -F "_" '{print $3}'); echo -e ${i}'\t'${sam}_${x}_${lane}_${read}_001.fastq.gz; done > ../names.txt  
mkdir ../files
cd ..
ln -s raw/* files/.
cd files
while read from to; do    echo "mv ${from} $to"; done < ../names.txt
```

## Steps

## 02 Read Trimming
sbatch 02_read_trimming.sh setup/02_read_trimming.sh

## 03 Mapping using BWA-mem
03_run_mapping.py -r [PATH/TO/REFERENCE] -i [PATH/TO/TRIMMED/FASTQ/FILES] -o [PATH/TO/OUTPUT/FOLDERS] -p [PATH/TO/THE/SCRIPTS/FOLDER]