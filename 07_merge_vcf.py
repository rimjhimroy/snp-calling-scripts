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
    parser.add_argument('-i', '--input', dest='infolder', help="full path of the folder with vcf files per sample [required]", required=True)
    parser.add_argument('-sc', '--scaff', dest='scaff', help="full path to the scaffold list to operate on [required]", required=True)
    parser.add_argument('-g', '--group', dest='group', help="genotyping group [required]", required=True)
    parser.add_argument('-s', '--script', dest='script', help="full path of scripts folder [required]", required=True)

    args = parser.parse_args()
    
    vcf=[f for f in os.listdir(args.infolder) if f.startswith(args.group) and f.endswith('.vcf.gz')]
    outname='%s_allchr.vcf.gz' % args.group 
    vcflist=["I="+f+" " for f in vcf]
    tomerge=' '.join(map(str, vcflist))
    cmd='picard-tools MergeVcfs %sO=%s' % (tomerge,outname)
    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)[1]