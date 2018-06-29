#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 5

import pandas as pd 
from subprocess import call
import sys
sra = pd.read_table('../data/SraRunTable.txt',header=None)
key1 = int(sys.argv[1])
key2 = int(sys.argv[2])

keys = sra.iloc[key1:key2,6]

#Fastq-dump commands from https://edwards.sdsu.edu/research/fastq-dump/
for key in keys:
	call(["fastq-dump",key,"gzip","--read-filter","pass","--outdir","../data/fastq","--skip-technical"])
