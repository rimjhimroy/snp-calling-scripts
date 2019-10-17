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
parser.add_argument('-i', '--input', dest='infolder', help="full path of the folder with input per library *_sort_mrkdup.fastq.gz files [required]", required=True)
parser.add_argument('-s', '--map', dest='map', help="full path to the file mapping library name to sample name [required]", required=True)
args = parser.parse_args()
