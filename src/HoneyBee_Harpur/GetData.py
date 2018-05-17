#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 5

import pandas as pd 
from subprocess import call
sra = pd.read_table('SraRunTable.txt')
keys = sra.iloc[:,6]

#Fastq-dump commands from https://edwards.sdsu.edu/research/fastq-dump/
for key in keys:
	call(["fastq-dump",key,"gzip","--read-filter","pass","--outdir","fastq","--skip-technical"])
