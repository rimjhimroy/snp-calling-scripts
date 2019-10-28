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
    parser.add_argument('-i', '--input', dest='infolder', help="full path of the folder with vcf files per sample of BQSR [required]", required=True)
    parser.add_argument('-r', '--ref', dest='ref', help="full path to the reference file [required]", required=True)
    parser.add_argument('-sc', '--scaff', dest='scaff', help="full path to the scaffold list to operate on [required]", required=True)
    parser.add_argument('-g', '--group', dest='group', help="genotyping group [required]", required=True)
    parser.add_argument('-s', '--script', dest='script', help="full path of scripts folder [required]", required=True)

    args = parser.parse_args()
    
    fa=args.ref
    if not os.path.isfile(fa):
            print("fasta file with reference genome is missing")
    refas=fa.rsplit('.', 1)[0]
    amb=args.ref+'.amb'
    ann=args.ref+'.ann'
    bwt=args.ref+'.bwt'
    pac=args.ref+'.pac'
    sa=args.ref+'.sa'
    fai=args.ref+'.fai'
    indexFiles=[amb, ann, bwt, pac, sa,fai]
    for i in indexFiles:
        if not os.path.isfile(i):
            print("Index files for reference genome are incomplete. Please check %s file" % i)
            break
        else:
            continue
    
    scriptpath=args.script
    basepath=args.infolder.rsplit('/', 1)[0]
    outfolder=basepath+"/joint_genotyping"
    mk_dir(outfolder)
    
    lncmd='ln -s %s %s/.' %(basepath+"/BQSR/*/*_final_variants.g.vcf.gz",outfolder+"/.")
    lnp = subprocess.Popen(lncmd, shell=True)
    lnsts = os.waitpid(lnp.pid, 0)[1]
    lncmd2='ln -s %s %s/.' %(basepath+"/BQSR/*/*_final_variants.g.vcf.gz.tbi",outfolder+"/.")
    lnp2 = subprocess.Popen(lncmd2, shell=True)
    lnsts2 = os.waitpid(lnp2.pid, 0)[1]
    
    vcf=[f for f in os.listdir(outfolder) if f.endswith('_final_variants.g.vcf.gz')]
    outname=set([i.split("_final_variants.g.vcf", 1)[0] for i in vcf])
    print(outname)
    GATK4="java -Xmx46g -Xms46g -Djava.oi.tmpdir=`pwd`/tmp -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar "
    
    loadmod="source %s" % basepath+"/modules.sh"
    vcflist=["-V "+f+"_final_variants.g.vcf.gz" for f in outname]
    tomerge=' '.join(map(str, vcflist))
    mk_dir(outfolder+"/scripts")
    with open(args.scaff) as fi:
        scaffs = fi.read().splitlines()
    for i in scaffs:
        scaff_name=i
        scaff_name=scaff_name.replace("|", "_")
        scaff_name=scaff_name.replace(".", "_")
        command1="%s GenomicsDBImport %s --reader-threads 8 --genomicsdb-workspace-path %s_%sdb -L %s" % (GATK4,tomerge,args.group,scaff_name,i.replace("|", "\|"))
        command2="%s GenotypeGVCFs --include-non-variant-sites -R %s -V gendb://%s_%sdb -O %s_%s.vcf.gz" % (GATK4,args.ref,args.group,scaff_name,args.group,scaff_name)

        filename='%s/joint_genotyping/scripts/joint_genotyping_%s.sh' %(basepath,scaff_name)
        script=open(filename, 'w')
        script.write('#!/bin/bash'+'\n')
        script.write(loadmod+'\n')
        script.write(command1+'\n')
        script.write(command2+'\n')
        script.close()
        
        cmd = ('sbatch -p all --mem=48G --time=5:00:00 '+filename)
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]