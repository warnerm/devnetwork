#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 5

import pandas as pd 
from subprocess import call
sra = pd.read_csv('drosophila_sra.csv')
keys = sra.iloc[:,0]

#Fastq-dump commands from https://edwards.sdsu.edu/research/fastq-dump/
for key in keys:
	call(["fastq-dump",key,"--split-3","--read-filter","pass","--outdir","fastq","--skip-technical"])
