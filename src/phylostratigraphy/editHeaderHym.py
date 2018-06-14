from Bio import SeqIO
import re, sys
import numpy as np

codes = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

prefix = infile.replace('.fa','')
prefix = prefix[-4:]
found = 0
with open(codes) as codes:
	for line in codes:
		if prefix in line:
			taxID = re.split('\t',line)[2].strip()
			found = 1

if found == 1:
	records = []
	seqiter = SeqIO.parse(open(infile), 'fasta') 
	for seq in seqiter:
		seq.description=taxID+seq.description
		records = np.append(records,seq)

	seqIO.write(records,outfile,"fasta")

else: #Just make an empty file if the genome isn't in the list
	out = open(outfile,'w')
	out.close()


