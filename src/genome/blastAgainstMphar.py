#!/usr/bin/python

#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 50

import sys,re
from joblib import Parallel, delayed
import numpy as np
import sqlite3
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
from subprocess import call

inputfile = sys.argv[1]
mphar = sys.argv[2]
output = mphar+"_blastmatch.txt"

seqiter = SeqIO.parse(open(inputfile),'fasta')

f = open(output,'w')
f.close()

def blastRecord(seq):
	num_match = 0

	SeqIO.write(seq, "temp"+seq.id+".fa", "fasta")

	tblastn_cline = NcbitblastnCommandline(query="temp"+seq.id+".fa",db=mphar,evalue=1e-10,
										outfmt=5,out="blast"+seq.id+".xml")
	stdout, stderr = tblastn_cline()

	result_handle = open("blast"+seq.id+".xml")
	blast_record = NCBIXML.read(result_handle)
	E_VALUE_THRESH = 1e-10	

	alignments = []

	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESH:
				alignments = np.append(alignments,str(alignment.title))

	alignment_unique = np.unique(alignments)

	with open(output,'a') as f:
		f.write(seq.id+'\t'+str(len(alignment_unique))+'\n')

	call(["rm","temp"+seq.id+".fa","blast"+seq.id+".xml"])


Parallel(n_jobs=50)(delayed(blastRecord)(seq) for seq in seqiter)

