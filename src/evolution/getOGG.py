#Get OGGs from a particular heirarchical category given a list of odb9 genes
from joblib import Parallel, delayed
import numpy as np, sys, re, pandas as pd


infile = sys.argv[1]
outfile = sys.argv[2]
threads = sys.argv[3]
OGGmap = sys.argv[4]
OGGlev = sys.argv[5]

def getMet(gene,outfile,OGGlev):
	with open(OGGmap) as ogMap:
		for line in ogMap:
			if gene in line:
				if OGGlev in line: #Restrict search to whatever level we care about
					with open(outfile,'a') as out:
						out.write(gene.strip()+'\t'+re.split('\t',line)[0]+'\n')

out = open(outfile,'w')
out.close()

with open(infile) as gene_file:
	Parallel(n_jobs=int(threads))(delayed(getMet)(gene_line,outfile,OGGlev) for gene_line in gene_file)
