#Get metazoan OGG given a list of odb9 genes
from joblib import Parallel, delayed
import numpy as np, sys, re, pandas as pd


infile = sys.argv[1]
outfile = sys.argv[2]
threads = sys.argv[3]
OGGmap = sys.argv[4]

def getMet(gene,outfile):
	#gene = re.split('\t',gene_line)[1].strip()
	with open(OGGmap) as ogMap:
		for line in ogMap:
			if gene in line:
				if "EOG090R" in line: #All metazoan OGGs start with this
					with open(outfile,'a') as out:
						out.write(gene.strip()+'\t'+re.split('\t',line)[0]+'\n')

out = open(outfile,'w')
out.close()

with open(infile) as gene_file:
	Parallel(n_jobs=int(threads))(delayed(getMet)(gene_line,outfile) for gene_line in gene_file)
