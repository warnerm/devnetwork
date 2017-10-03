#!/usr/local/bin/python2.7

from Bio import SeqIO
import sys, getopt
from subprocess import call
from joblib import Parallel, delayed

def InOut(argv):
	inputfile=''
	outputfile = ''
	seqfile = ''
	try:
	        opts, args = getopt.getopt(argv,"hi:o:s:",["ifile=","ofile=","sfile="])
	except getopt.GetoptError:
	        print 'KeepSeq.py -i <inputfile> -o <outputfile> -s <seqfile> -t <threads>'
	        sys.exit(2)
	for opt, arg in opts:
	        if opt in ("-i","--ifile"):
	                inputfile = arg
	        elif opt in ("-o","--ofile"):
	                outputfile = arg   
	        elif opt in ("-s","--ofile"):
	                seqfile = arg                  
	return inputfile,outputfile,seqfile

def main(argv):
	inputfile,outputfile,seqfile = InOut(argv)
	wanted = [line.strip() for line in open(seqfile)]                               
	seqiter = SeqIO.parse(open(inputfile), 'fasta')                                    
	SeqIO.write((seq for seq in seqiter if seq.description in wanted), outputfile, "fasta")
		    
if __name__ == "__main__":
	main(sys.argv[1:]) 

