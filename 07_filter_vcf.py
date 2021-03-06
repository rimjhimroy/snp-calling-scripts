#! /usr/bin/env python3
# Created by Rimjhim Roy Choudhury

import argparse, time, os, subprocess, errno
import pandas as pd


def mk_dir(path):
    try:
        os.mkdir(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

if __name__ == '__main__':
    # create variables that can be entered as arguments in command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', help="full path of the merged  and genotyped vcf file to filter [required]", required=True)
    parser.add_argument('-r', '--reference', dest='r', help="full path of the ref genome index files [required]", required=True)

    args = parser.parse_args()
    GATK4="java -Xmx46g -Xms46g -Djava.oi.tmpdir=`pwd`/tmp -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar "
    
    basepath=args.input.rsplit('/', 1)[0]
    loadmod="source %s" % basepath+"/modules.sh"
    
    baseout=os.path.splitext(os.path.basename(args.input))[0]
    outname1=baseout+"_raw_SNPS.vcf.gz"
    outname2=baseout+"_filtered_SNPS.vcf.gz"
    cmd1 = '%s SelectVariants -R %s -V %s -select-type SNP -select-type NO_VARIATION -O %s' % (GATK4,args.r,args.input,outname1)
    cmd2 = "%s VariantFiltration -R %s -V %s -O %s -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'"  % (GATK4,args.r,outname1,outname2)
 
    outname3=baseout+"_raw_INDEL.vcf.gz"
    outname4=baseout+"_filtered_INDEL.vcf.gz"
    cmd3 = '%s SelectVariants -R %s -V %s -select-type INDEL -O %s' % (GATK4,args.r,args.input,outname3)
    cmd4 = "%s VariantFiltration -R %s -V %s -O %s -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20'" % (GATK4,args.r,outname3,outname4)
 
    filename1='%s/filtering_SNP.sh' %(basepath)
    script=open(filename1, 'w')
    script.write('#!/bin/bash'+'\n')
    script.write(loadmod+'\n')
    script.write(cmd1+'\n')
    script.write(cmd2+'\n')
    script.close()
    
    scmd = ('sbatch -p all --mem=48G --time=5:00:00 '+filename1)
    sp = subprocess.Popen(scmd, shell=True)
    ssts = os.waitpid(sp.pid, 0)[1]
      
    filename2='%s/filtering_INDELs.sh' %(basepath)  
    script=open(filename2, 'w')
    script.write('#!/bin/bash'+'\n')
    script.write(loadmod+'\n')
    script.write(cmd3+'\n')
    script.write(cmd4+'\n')
    script.close()
    
    icmd = ('sbatch -p all --mem=48G --time=5:00:00 '+filename2)
    ip = subprocess.Popen(icmd, shell=True)
    ists = os.waitpid(ip.pid, 0)[1]