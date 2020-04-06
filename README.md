# Snp calling scripts

## Preprocess

You should not merge the different FASTQs from different lanes before aligning them. If you do so, then you will loose "read group" (RG tag) information. I would align each lane separately, then the scripts add specific RG tags to each lane and then merge the aligned bams (respecting RG information). RG information is important for downstream analyses in GATK, such as  BQSR.  

To get help about how to run the python scripts launch: `scriptname.py -h`

1. Rename [SAMPLENAME]\_[LANE]\_[PAIR]\_[random_string].fastq.gz files by adding the run id to the file name, giving [SAMPLE]\_[LANE]\_[PAIR]\_001.fastq.gz (important if there are several runs from the same library, else ignore). Make sure there is no '_' in SAMPLENAME.  

```bash
cd raw # go to your raw directory
for i in $( ls *.gz ); do x=$(zcat $i | head -n 1 | awk -F ":" '{print $2}');  sam=$(echo $i | awk -F "_" '{print $1}'); lane=$(echo $i | awk -F "_" '{print $2}');read=$(echo $i | awk -F "_" '{print $3}'); echo -e ${i}'\t'${sam}_${x}_${lane}_${read}_001.fastq.gz; done > ../names.txt  #create a table to rename files
mkdir ../files
cd ../files
ln -s ../raw/* .  #create softlinks to the files
while read from to; do mv ${from} $to; done < ../names.txt #rename the files
```

2. FastQC and multiqc

'setup/run_qc.json' is a setup file that you need to change with the location of your data before running  
It will run fastqc on the fastq files and generate quality plots which is useful to spot potential problems in high througput sequencing datasets. It highlight any areas where this library looks unusual (highlighting them in red) and where you should take a closer look.  

```bash
sbatch doQC.sh setup/run_qc.json
```

## Steps

## 01 Read Trimming

'setup/01_read_trimming.json' is a setup file that you need to change with the location of your data before running  
When data is sequenced on Illumina, adapters are added for the fragments to attach to the beads. If these adapters are not removed they can result in false assembly or other issues. Trimmomatic is used in the script to remove sequence adapters.  

```bash
sbatch 01_read_trimming.sh setup/01_read_trimming.json
```

## 02 Mapping using BWA-mem

```bash
usage: 02_run_mapping.py [-h] -r R -i INFOLDER -o OUTPUT -p SCRIPT
                         [-N THREADS] [--mem MEM] [--time TIME] [-k SEED]
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
  -v VERSION, --version VERSION
                        bwa version to use [default 0.7.17]
  --samtools SAMTOOLS   samtools version [default 1.8]
  --print PRT           If changed to true then shell files are printed to
                        screen and not launched [false]
```

## 03 Merge different bam files from different libraries for the same sample

```bash
usage: 03_merge_libraries_into_sample.py [-h] -i INFOLDER -s MAP

optional arguments:
  -h, --help            show this help message and exit
  -i INFOLDER, --input INFOLDER
                        full path of the folder with input per library
                        *_sort_mrkdup.bam files [required]
  -s MAP, --map MAP     full path to the file mapping library name to sample
                        name [required]
```

## 04 Run BQSR and call variants

```bash
usage: 04_run_variant_calling.py [-h] -i INFOLDER -r REF -s SCRIPT -p PLOIDY

optional arguments:
  -h, --help            show this help message and exit
  -i INFOLDER, --input INFOLDER
                        full path of the folder with merged bam files
                        [required]
  -r REF, --ref REF     full path to the reference file [required]
  -s SCRIPT, --script SCRIPT
                        full path of scripts folder [required]
  -p PLOIDY, --ploidy PLOIDY
                        full path of the file containing sample to ploidy map
                        [required]
```

## 05 Joint genotyping

```bash
usage: 05_joint_genotype_vcf.py [-h] -i INFOLDER -r REF -sc SCAFF -g GROUP -s
                                SCRIPT

optional arguments:
  -h, --help            show this help message and exit
  -i INFOLDER, --input INFOLDER
                        full path of the folder with vcf files per sample of
                        BQSR [required]
  -r REF, --ref REF     full path to the reference file [required]
  -sc SCAFF, --scaff SCAFF
                        full path to the scaffold list to operate on
                        [required]
  -g GROUP, --group GROUP
                        genotyping group [required]
  -s SCRIPT, --script SCRIPT
                        full path of scripts folder [required]
```

## 06 Merge VCFs from different chromosomes and scaffolds

```bash
usage: 06_merge_vcf.py [-h] -i INFOLDER -sc SCAFF -g GROUP

optional arguments:
  -h, --help            show this help message and exit
  -i INFOLDER, --input INFOLDER
                        full path of the folder with vcf files per chromosome
                        [required]
  -sc SCAFF, --scaff SCAFF
                        full path to the scaffold list to operate on
                        [required]
  -g GROUP, --group GROUP
                        genotyping group [required]
```

## 07_filter_vcf.py

```bash
usage: 07_filter_vcf.py [-h] -i INPUT -r R

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        full path of the merged and genotyped vcf file to
                        filter [required]
  -r R, --reference R   full path of the ref genome index files [required]
```
