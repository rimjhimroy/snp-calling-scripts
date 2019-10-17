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

class sampmap:
    def __init__(self,entry):
        popmap=pd.read_csv(entry, delimiter='\s+',lineterminator='\n',chunksize =50)
        map_df = pd.concat(popmap, ignore_index=True)
        map_df.columns = ['lane', 'sample']
        self.lane = list(map_df['lane'])
        sample = list(map_df['sample'])
        self.map = dict(zip(self.lane, sample))
        self.sample = set(sample)



if __name__ == '__main__':
    # create variables that can be entered as arguments in command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='infolder', help="full path of the folder with input per library *_sort_mrkdup.fastq.gz files [required]", required=True)
    parser.add_argument('-s', '--map', dest='map', help="full path to the file mapping library name to sample name [required]", required=True)
    args = parser.parse_args()

    map_values=sampmap(args.map)
    lanedict=map_values.map
    samples=map_values.sample
    basedir=args.infolder.rsplit('/', 1)[0]
    mk_dir(basedir+"/mergedbam")
    for i in samples:
        lanelist=[basedir+"/mapout/"+lane+"_sort_mkrdup.bam" for lane, slist in lanedict.items() if slist==i]
        tomerge=' '.join(map(str, lanelist))
        outbam=i+".bam"

        cmd='sbatch -p all --wrap "samtools merge %s/mergedbam/%s %s"' %(basedir,outbam,tomerge)
        print(cmd)
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]