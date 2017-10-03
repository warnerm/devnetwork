#OGGcreate.py create OGG files; python 3.5

import numpy as np, sys, getopt, pandas as pd
from subprocess import call
from joblib import Parallel, delayed
import re

def InOut(argv):
	inputfile=''
	outputfile = ''
	try:
	        opts, args = getopt.getopt(argv,"hi:t:s:")
	except getopt.GetoptError:
	        print ('OGGcreate.py -i <inputfile> -s <seqfile>')
	        sys.exit(2)
	for opt, arg in opts:
	        if opt in ("-i","--ifile"):
	                inputfile = arg  
	        elif opt in ("-t","--tfile"):
	                num_cores = arg  
	        elif opt in ("-s","--sfile"):
	                seq_file = arg                     
	return inputfile,num_cores,seq_file

def createFasta(OGG,data,inputfile):
Odat = data[data[:,1]==OGG,:] #filter ODB8 file for the specific OGG
outputfile=OGG+".fa"
get = ''
for row in range(np.shape(Odat)[0]):
	line = '\t'.join(str(e) for e in Odat[row,])
	gene = re.sub(r'.*\t([A-Z0-9]+)\t([A-Za-z]*).*\t([0-9]+:[0-9a-zA-Z]+).*',r'\g<1>|\g<2>|\g<3>',line)
	get = get + ' ' + gene
call(["pyfasta","extract","--header","--fasta",inputfile,"'EOG8004CG|Lalb|88501:0018c5'"],stderr=outputfile)

def main(argv):
inputfile,num_cores,seq_file = InOut(argv)
df = pd.read_csv(seq_file,sep="\t")
data = np.array(df)
OGGS = np.unique(data[:,1])
	Parallel(n_jobs=int(num_cores))(delayed(createFasta)(OGG,data,inputfile)for OGG in OGGS)
