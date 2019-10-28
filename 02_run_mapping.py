#! /usr/bin/env python3
# Created by Rimjhim Roy Choudhury

import argparse, time, os, subprocess, errno

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
parser.add_argument('-i', '--input', dest='infolder', help="full path of the folder with input *.trim.fastq.gz files [required]", required=True)
parser.add_argument('-o', '--output', dest='output', help="full path of output folder [required]", required=True)
parser.add_argument('-p', '--script', dest='script', help="full path of scripts folder [required]", required=True)
parser.add_argument('-N', '--threads', dest='threads', help="number of threads [default 4]", default=4, type=int)
parser.add_argument('--mem', dest='mem', help="memory [default 16G]", default="16G")
parser.add_argument('--time', dest='time', help="maximum run time as hours:minutes:seconds [default 4:0:0]", default="24:0:0")
parser.add_argument('-k', '--seed', dest='seed', help="minimum seed length [default 19]", default=19, type=int)
parser.add_argument('-v', '--version', dest='version', help="bwa version to use [default 0.7.17]", default="0.7.17")
parser.add_argument('--samtools', dest='samtools', help="samtools version [default 1.8]", default='1.8')
parser.add_argument('--print', type=str, dest='prt', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')


args = parser.parse_args()

#scriptpath=script.rsplit('/', 1)[0]
scriptpath=args.script
fastq=[f for f in os.listdir(args.infolder) if f.endswith('.trim.fastq.gz')]
sfq=[i for i in fastq if "_L1_" in i] 
outname=set([i.split("_L1_", 1)[0] for i in sfq])
d = dict(zip(outname,outname))
for key, value in d.items():
    d[key] = value.split('_')[0]
mk_dir(args.output+"/mapout")
mk_dir(args.output+"/maplog")
mk_dir(args.output+"/mapscripts")
for key in d:
    pair1lane1=args.infolder+"/"+key+"_L1_R1_001.trim.fastq.gz"
    pair2lane1=args.infolder+"/"+key+"_L1_R2_001.trim.fastq.gz"
    pair1lane2=args.infolder+"/"+key+"_L2_R1_001.trim.fastq.gz"
    pair2lane2=args.infolder+"/"+key+"_L2_R2_001.trim.fastq.gz"
    outfileL1=args.output+"/mapout/"+key+"_L1"
    outfileL2=args.output+"/mapout/"+key+"_L2"
    jobnameL1=key+"_L1"
    jobnameL2=key+"_L2"
    cmdL1='%s/BWA-mem.py -r %s -i %s -m %s --sample %s --library %s -o %s -n %s -N %s --mem %s --time %s -k %s -v %s --samtools %s --print %s' %(scriptpath,args.r,pair1lane1,pair2lane1,d[key],d[key],outfileL1,jobnameL1,args.threads,args.mem,args.time,args.seed,args.version,args.samtools,args.prt)
    cmdL2='%s/BWA-mem.py -r %s -i %s -m %s --sample %s --library %s -o %s -n %s -N %s --mem %s --time %s -k %s -v %s --samtools %s --print %s' %(scriptpath,args.r,pair1lane2,pair2lane2,d[key],d[key],outfileL2,jobnameL2,args.threads,args.mem,args.time,args.seed,args.version,args.samtools,args.prt) 
    print(cmdL1)
    pL1 = subprocess.Popen(cmdL1, shell=True)
    stsL1 = os.waitpid(pL1.pid, 0)[1]
    print(cmdL2)
    pL2 = subprocess.Popen(cmdL2, shell=True)
    stsL2 = os.waitpid(pL2.pid, 0)[1]



