#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 20

import numpy as np, sys, re, getopt, pandas as pd
import subprocess
from subprocess import call
from joblib import Parallel, delayed
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML

def InOut(argv):
	try:
	        opts, args = getopt.getopt(argv,"h:s:d:t:o:",["sfile=","dfile=","tfile=","ofile="])
	except getopt.GetoptError:
	        print 'pairwiseBlast.py -s <seqfile> -d <db> -t <threads>'
	        sys.exit(2)
	for opt, arg in opts:
	        if opt in ("-s","--sfile"):
	                seqfile = arg
	        if opt in ("-d","--dfile"):
	                dbfile = arg 
	        if opt in ("-t","--tfile"):
	                threads = arg   
	        if opt in ("-o","--ofile"):
	                outfile = arg                         
	return seqfile,dbfile,threads,outfile

def callProc(gene,seqfile,dbfile,outfile):
	spec = seqfile.replace('_prot.fa','')

	#Make fasta file of the individual protein
	seqiter = SeqIO.parse(open(seqfile),'fasta')
	SeqIO.write((seq for seq in seqiter if seq.id in gene), "temp"+gene+".fa", "fasta")

	blastp_cline = NcbiblastpCommandline(query="temp"+gene+".fa",db=dbfile,evalue=1e-10,
											outfmt=5,out="blast"+gene+".xml")

	stdout, stderr = blastp_cline()

	result_handle = open("blast"+gene+".xml")
	blast_record = NCBIXML.read(result_handle)
	E_VALUE_THRESH = 1e-10

	#Get first match
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESH:
				with open(outfile,'a') as f:
					f.write(gene+'\t'+alignment.title+'\n')
				break
		break

	#Remove temporary files
	call(["rm","temp"+gene+".fa","blast"+gene+".xml"])

def main(argv):
	seqfile,dbfile,threads,outfile = InOut(argv)
	dbfile = dbfile.replace('.phr','')
	seq_IDs = []
	out = open(outfile,'w')
	out.close()

	#Get list of all proteins
	with open(seqfile,'r') as f:
		for line in f:
			if ">" in line:
				l = line.replace(">","")
				l = l.replace("\n","")
				seq_IDs = np.append(seq_IDs,l)

	Parallel(n_jobs=int(threads))(delayed(callProc)(seq,seqfile,dbfile,outfile) for seq in seq_IDs)

if __name__ == "__main__":
	main(sys.argv[1:])