import sys, getopt
from Bio.Align.Applications import ClustalwCommandline
import re, os
import numpy as np
from Bio.Phylo.PAML import codeml
from subprocess import call
from joblib import Parallel, delayed
from Bio import SeqIO

import shutil

global seq1,seq2,seq1_nucl,seq2_nucl
mapfile = sys.argv[1]
seq1 = sys.argv[2]
seq2 = sys.argv[3]
seq1_nucl = sys.argv[4]
seq2_nucl = sys.argv[5]
outfile = sys.argv[6]
threads = sys.argv[7]

with open(outfile,'w') as out:
	out.write('Gene\tdN_dS\n')

def alignGene(line):
	try:
		gene1 = line.split(" ")[2]
		gene2 = line.split(" ")[7]
		print(gene1)
		print(gene2)
		cds1 = line.split(" ")[1]
		cds2 = line.split(" ")[6]
		geneName = line.split(" ")[4]
		#Make file with both protein sequences
		seqiter = SeqIO.parse(open(seq1),'fasta')
		SeqIO.write((seq for seq in seqiter if seq.id == gene1), "scratch/"+gene1+".fa", "fasta")
		seqiter = SeqIO.parse(open(seq2),'fasta')
		SeqIO.write((seq for seq in seqiter if seq.id == gene2), "scratch/"+gene2+".fa", "fasta")

		with open("scratch/"+geneName+".fa",'wb') as wfd:
			for f in ["scratch/"+gene1+".fa","scratch/"+gene2+".fa"]:
				with open(f,'rb') as fd:
					shutil.copyfileobj(fd, wfd) 

		#Make file with both protein sequences
		seqiter = SeqIO.parse(open(seq1_nucl),'fasta')
		SeqIO.write((seq for seq in seqiter if seq.id == cds1), "scratch/"+cds1+".fa", "fasta")
		seqiter = SeqIO.parse(open(seq2_nucl),'fasta')
		SeqIO.write((seq for seq in seqiter if seq.id == cds2), "scratch/"+cds2+".fa", "fasta")

		with open("scratch/"+geneName+"_nucl.fa",'wb') as wfd:
			for f in ["scratch/"+cds1+".fa","scratch/"+cds2+".fa"]:
				with open(f,'rb') as fd:
					shutil.copyfileobj(fd, wfd) 

		#Make tree file
		cline = ClustalwCommandline("clustalw2", 
			infile="scratch/"+geneName+".fa",
			newtree='scratch/'+geneName+'tree.tre')
		stdout, stderr = cline()

		#Make alignment file
		cline = ClustalwCommandline("clustalw2", 
			infile="scratch/"+geneName+".fa",
			output="CLUSTAL",
			outfile='scratch/'+geneName+'alignment.aln')
		stdout, stderr = cline()

		cmd = "perl pal2nal.pl scratch/"+geneName+"alignment.aln scratch/"+geneName+"_nucl.fa -output paml > scratch/"+geneName+"alignment_nucl.phy"
		#Run pal2nal
		os.system(cmd) 
		cml = codeml.Codeml()
		cml.alignment = 'scratch/'+geneName+'alignment_nucl.phy'
		cml.working_dir = "./scratch"
		cml.tree = 'scratch/'+geneName+'tree.tre'
		cml.out_file = 'scratch/'+geneName+'out.txt'
		cml.set_options(seqtype=1,
	   	        verbose=1,
	            noisy=0,
	            model=1,
	            runmode=-2,
	            Mgene=0,
	            NSsites=[0],
	            CodonFreq=2,
	            cleandata=1)
		cml.run(verbose=False)
		
		with open(cml.out_file) as results:
			with open(outfile,'a') as out:
				for line in results:
					if line.find("dN/dS=") > -1:
						line = line.split()
						out.write(geneName+'\t'+str(line[line.index("dN/dS=")+1])+'\n')
		cmd = "rm scratch/*"+geneName+"* scratch/"+cds1+"* scratch/"+cds2+"* scratch/*"+gene1+"* scratch/*"+gene2+"*"
		os.system(cmd)
	except:
		return



with open(mapfile,'r') as m:
	Parallel(n_jobs=int(threads))(delayed(alignGene)(line) for line in m)

