#! /usr/bin/env python3
# Created by Rimjhim Roy Choudhury

import argparse, time, os, subprocess, errno
import pandas as pd

class ploidymap:
    def __init__(self,entry):
        popmap=pd.read_csv(entry, delimiter='\s+',lineterminator='\n',chunksize =50,header=None)
        map_df = pd.concat(popmap, ignore_index=True)
        map_df.columns = ['sample', 'ploidy']
        self.sample = list(map_df['sample'])
        ploidy = list(map_df['ploidy'])
        self.map = dict(zip(self.sample, ploidy))
        self.ploidy = set(ploidy)

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
    parser.add_argument('-i', '--input', dest='infolder', help="full path of the folder with merged bam files [required]", required=True)
    parser.add_argument('-r', '--ref', dest='ref', help="full path to the reference file [required]", required=True)
    parser.add_argument('-s', '--script', dest='script', help="full path of scripts folder [required]", required=True)
    parser.add_argument('-p', '--ploidy', dest='ploidy', help="full path of the file containing sample to ploidy map [required]", required=True)

    args = parser.parse_args()
    basedir=args.infolder.rsplit('/', 1)[0]
    mk_dir(basedir+"/BQSR")
    mk_dir(basedir+"/BQSRlog")
    logoutput='%s/BQSRlog/BQSR-%%N-%%j.log' %(basedir)
    logerror='%s/BQSRlog/BQSR-%%N-%%j.log' %(basedir)
    bam=[f for f in os.listdir(args.infolder) if f.endswith('.bam')]
    map_values=ploidymap(args.ploidy)
    pdict=map_values.map
    samples=map_values.sample
    
    for i in bam:
        sam=i.rsplit('.', 1)[0]
        ploidy=pdict[sam]
        sample=basedir+"/mergedbam/"+i
        cmd='sbatch -p all -J %s --output=%s --error=%s %s/variant_calling.sh %s %s %s %s' %(sam,logoutput,logerror,args.script,args.ref,sample,basedir,ploidy)
        print(cmd)
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
        print(cmd)