from Bio import SeqIO
import sys
# take the output of gffread and find isoform length

genes = SeqIO.to_dict(SeqIO.parse(sys.argv[1],"fasta"))

with open(outfile,'a') as out:
	out.write("GeneID\tLength\n")
	for rec in genes:
		out.write(genes[rec].description.split("=")[0]+'\t'+str(len(rec))+'\n')
	