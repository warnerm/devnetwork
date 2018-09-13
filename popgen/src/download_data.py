from subprocess import call
import sys
sra = sys.argv[1]

call(["fastq-dump",sra,"gzip","--split-3","--read-filter","pass","--outdir","../data/fastq","--skip-technical"])
