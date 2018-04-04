#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 1

import pandas as pd 
from subprocess import call
sra = pd.read_table('samples_jasper.txt')
keys = sra.iloc[:,3]

#Fastq-dump commands from https://edwards.sdsu.edu/research/fastq-dump/
for key in keys:
	call(["fastq-dump",key,"gzip","--split-3","--read-filter","pass","--outdir","fastq","--skip-technical"])
