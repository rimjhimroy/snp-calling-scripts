# Snp calling scripts

## Preprocess

You should not merge the different FASTQs from different lanes before aligning them. If you do so, then you will loose "read group" (RG tag) information. I would align each lane separately, then the scripts add specific RG tags to each lane and then merge the aligned bams (respecting RG information). RG information is important for downstream analyses in GATK, such as  BQSR.  

To get help about how to run the python scripts launch: `scriptname.py -h`

1. Rename [SAMPLENAME]\_[LANE]\_[PAIR]\_[random_string].fastq.gz files by adding the run id to the file name, giving [SAMPLE]\_[LANE]\_[PAIR]\_001.fastq.gz (important if there are several runs from the same library, else ignore). Make sure there is no '_' in SAMPLENAME.  

```bash
cd raw
for i in $( ls *.gz ); do x=$(zcat $i | head -n 1 | awk -F ":" '{print $2}');  sam=$(echo $i | awk -F "_" '{print $1}'); lane=$(echo $i | awk -F "_" '{print $2}');read=$(echo $i | awk -F "_" '{print $3}'); echo -e ${i}'\t'${sam}_${x}_${lane}_${read}_001.fastq.gz; done > ../names.txt  
mkdir ../files
cd ..
ln -s raw/* files/.
cd files
while read from to; do    echo "mv ${from} $to"; done < ../names.txt
```

2. FastQC and multiqc

'setup/run_qc.json' is a setup file that you need to change with the location of your data before running

```bash
sbatch doQC.sh setup/run_qc.json
```

## Steps

## 01 Read Trimming

'setup/01_read_trimming.sh' is a setup file that you need to change with the location of your data before running

```bash
sbatch 01_read_trimming.sh setup/01_read_trimming.sh
```

## 02 Mapping using BWA-mem

```bash
02_run_mapping.py -r [PATH/TO/REFERENCE] -i [PATH/TO/TRIMMED/FASTQ/FILES] -o [PATH/TO/OUTPUT/FOLDERS] -p [PATH/TO/THE/SCRIPTS/FOLDER]
```

```bash
usage: 02_run_mapping.py [-h] -r R -i INFOLDER -o OUTPUT -p SCRIPT
                         [-N THREADS] [--mem MEM] [--time TIME] [-k SEED] [-s]
                         [-v VERSION] [--samtools SAMTOOLS] [--print PRT]

optional arguments:
  -h, --help            show this help message and exit
  -r R, --reference R   full path of the ref genome index files [required]
  -i INFOLDER, --input INFOLDER
                        full path of the folder with input *.trim.fastq.gz
                        files [required]
  -o OUTPUT, --output OUTPUT
                        full path of output folder [required]
  -p SCRIPT, --script SCRIPT
                        full path of scripts folder [required]
  -N THREADS, --threads THREADS
                        number of threads [default 4]
  --mem MEM             memory [default 16G]
  --time TIME           maximum run time as hours:minutes:seconds [default
                        4:0:0]
  -k SEED, --seed SEED  minimum seed length [default 19]
  -s, --no_scratch      don\'t use local scratch. CURRENTLY CANNOT USE SCRATCH
                        [default DON\'T use local scratch]
  -v VERSION, --version VERSION
                        bwa version to use [default 0.7.17]
  --samtools SAMTOOLS   samtools version [default 1.8]
```
## 03 Merge different bam files from different libraries for the same sample
