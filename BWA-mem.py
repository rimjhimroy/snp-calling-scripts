#! /usr/bin/env python3

# script uses bwa mem which is recommended for Illumina reads >70bp


## CANNOT use scratch with slurm at the moment

#### Necessary input ####
# the script assumes that index files for the reference genome have been created before (bwa index -a bwtsw <ref.fa>)
# and that they are in the same folder as the fasta file with the reference genome
# for paired-end reads, first and second read are in different files (but at same position in the file)

#### Mapping parameters ####
# seed length is currently the only parameter that can be set. All other mapping parameters are default
# bwa mem performs local alignment which means that there can be multiple primary alignments for a read if different parts of the read
# map to different places. -M means that only one of the alignments is indicated as primary, all others as secondary. This is required
# by Picard tools

#### Notes on bwa mem ####
# The manual and paper (Li 2013 unpublished) are really unclear. If I understand correctly, bwa estimates insert length from the data
# it then decides whether two reads should be output as paired or unpaired depending on a score that is influenced by a penalty
# given for reporting unpaired read pairs and a penalty for having a 'strange' insert size.

# bwa reports end-to-end alignment if this is better than the best score for a local alignment - a "local alignment penalty" (clipping penalty -L)

# Base quality is not considered

### output ###
# currently only the best hit(s) are output (I think). To output all found alignments (above a certain quality threshold?), set -a



import argparse, time, os, gzip, re, subprocess, errno

def mk_dir(path):
    try:
        os.mkdir(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)


parser = argparse.ArgumentParser()

parser.add_argument('-r', '--reference', dest='r', help="full path of the ref genome index files [required]", required=True)
parser.add_argument('-i', '--input', dest='mate1', help="full path fastq(.gz) file with first read [required]", required=True)
parser.add_argument('-m', '--mate', dest='mate2', help="full path fastq(.gz) file with second read [optional]")
parser.add_argument('-o', '--output', dest='output', help="full path of output file [required]", required=True)
parser.add_argument('--sample', dest='sample', help="sample name to be written to readgroup field SM [required]", required=True)
parser.add_argument('--library', dest='library', help="library name to be written to readgroup field LB [require]", required=True)
parser.add_argument('-N', '--threads', dest='threads', help="number of threads [default 4]", default=4, type=int)
parser.add_argument('--mem', dest='mem', help="memory [default 16G]", default="16G")
parser.add_argument('--time', dest='time', help="maximum run time as hours:minutes:seconds [default 4:0:0]", default="4:0:0")
parser.add_argument('-k', '--seed', dest='seed', help="minimum seed length [default 19]", default=19, type=int)
parser.add_argument('-n', '--name', dest='n', help="job name [default bwa]", default="bwa")
parser.add_argument('-s', '--no_scratch', dest='s', action='store_false', default=False, help="don't use local scratch. CURRENTLY CANNOT USE SCRATCH [default DON'T use local scratch]")
parser.add_argument('-v', '--version', dest='version', help="bwa version to use [default 0.7.17]", default="0.7.17")
parser.add_argument('--samtools', dest='samtools', help="samtools version [default 1.8]", default='1.8')
parser.add_argument('--print', type=str, dest='prt', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')


args = parser.parse_args()

# perform some checks to see if provided arguments are valid
fa=args.r
if not os.path.isfile(fa):
        print("fasta file with reference genome is missing")
refas=fa.rsplit('.', 1)[0]
amb=args.r+'.amb'
ann=args.r+'.ann'
bwt=args.r+'.bwt'
pac=args.r+'.pac'
sa=args.r+'.sa'
indexFiles=[amb, ann, bwt, pac, sa]
for i in indexFiles:
        if not os.path.isfile(i):
                print("Index files for reference genome are incomplete. Please check")
                break
        else:
                continue
if not os.path.isfile(args.mate1):
        print("fastq file with first read is missing")
if args.mate2:
        if not os.path.isfile(args.mate2):
                print("fastq file with second read is missing")



# Check if the fastq file for mate1 is zipped or not
if re.search('.gz$', args.mate1):
        zipped='yes'
else:
        zipped='no'

#set up read group
if zipped=='yes':
        mate1=gzip.open(args.mate1, 'r')
else:
        zipped='no'

#set up read group
if zipped=='yes':
        mate1=gzip.open(args.mate1, 'r')
else:
        mate1=open(args.mate1, 'r')
header=str(mate1.readline())
mate1.close()
header=header.split(":")
flowcell=header[2]
lane=header[3]
readgroupID=[flowcell,lane]
readgroupID='-'.join(readgroupID)

basefolder=args.output.rsplit('/', 2)[0]

# set up shell script
filename='%s/mapscripts/BWA_%s.sh' %(basefolder,args.n)
script=open(filename, 'w')
script.write('#!/bin/bash'+'\n')
script.write('#SBATCH --job-name=%s' %args.n +'\n')
script.write('#SBATCH --ntasks=1' +'\n')
script.write('#SBATCH --cpus-per-task=%i' %args.threads + '\n')
script.write('#SBATCH --output=%s/maplog/%s.out' %(basefolder,args.n) + '\n')
script.write('#SBATCH --error=%s/maplog/%s.err' %(basefolder,args.n) + '\n')
script.write('#SBATCH --mem=%s' %args.mem + '\n')
script.write('#SBATCH --time=%s' %args.time + '\n')


script.write('/bin/echo Job runs on host `hostname`'+'\n')
script.write('/bin/echo Job start time: `date`'+'\n')

script.write('export PATH=/software/bin:$PATH'+'\n')
script.write('module use /software/module/'+'\n')
script.write('module add vital-it'+'\n')
script.write('module add UHTS/Aligner/bwa/%s' %args.version +'\n')
script.write('module add UHTS/Analysis/samtools/%s' %args.samtools + '\n')
script.write('module add UHTS/Analysis/picard-tools/2.18.11'+'\n')

# change to local scratch (unless -s is provided)
if args.s:
        script.write('home=`pwd`'+'\n')    # absolut current path

#copy all required files to scratch:
        script.write('cp %s.* /scratch/"$JOB_ID".1.all.q/' %args.r +'\n') #copy ref genome and all index files
        script.write('cp %s /scratch/"$JOB_ID".1.all.q/' %args.mate1 +'\n') #copy first fastq file
        if args.mate2:
                script.write('cp %s /scratch/"$JOB_ID".1.all.q/' %args.mate2 +'\n') #copy second fastq file

        script.write('cd /scratch/"$JOB_ID".1.all.q/' +'\n')
        script.write('work=job_%s' %args.n + '\n')
        script.write('mkdir -p $work' + '\n')
        script.write('cd $work' + '\n')

# run bwa
script.write('bwa mem ' )
script.write('-t %i ' %args.threads)
script.write('-k %i ' %args.seed)
script.write("-R '@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\\tLB:%s\\tPU:%s-%s' " %(readgroupID, args.sample, args.library,readgroupID,args.sample)) #readgroups
script.write('-M ') # output only one alignment as primary (required for picard tools)
if args.s:
        refname=args.r
        refname=refname.split('/')
        refname=refname[len(refname)-1]
        script.write('../%s ' %refname)
        mate1=args.mate1
        mate1=mate1.split('/')
        mate1=mate1[len(mate1)-1]
        script.write('../%s ' %mate1)
        if args.mate2:
                mate2=args.mate2
                mate2=mate2.split('/')
                mate2=mate2[len(mate2)-1]
                script.write('../%s ' %mate2)
else:
        script.write('%s ' %args.r)
        script.write('%s ' %args.mate1)
        if args.mate2:
                script.write('%s ' %args.mate2)

script.write('samtools view -1 - > %s.bam' %args.output +'\n')



# bring back output if local scratch was used
if args.s:
        script.write('cd ..' + '\n')
        script.write('mv $work $home' + '\n')
        script.write('rm -rf $work' + '\n')
        script.write('cd $home' + '\n')

#Sort and mark duplicates with picard tools
outbam=args.output.rsplit('.', 1)[0]+"_sort.bam"
mrkdup=args.output.rsplit('.', 1)[0]+"_sort_mrkdup.bam"
tmpdir='%s/tmp/%s/' %(basefolder,args.n)
mk_dir(basefolder+"/tmp")
mk_dir(tmpdir)
#mk_dir(basefolder+"/mrkduplog")
metrics=basefolder+"/mrkduplog/"+args.n+"_metrics.txt"
script.write("picard-tools SortSam INPUT=%s OUTPUT=%s SORT_ORDER=coordinate TMP_DIR=%s" %(args.output+".bam", outbam,tmpdir)+'\n')
script.write("rm  %s" %(args.output+".sam")+'\n')
script.write("picard-tools MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s TMP_DIR=%s" %(outbam, mrkdup,metrics,tmpdir)+'\n')

script.write('/bin/echo Job end time: `date`'+'\n')
script.close()

if args.prt=='false':
        #send slurm job to NBI SLURM cluster
        cmd = ('sbatch '+filename)
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
else:
        file = open(filename,'r')
        data = file.read()
        print(data)

