from Bio import SeqIO
import sys
# take the output of gffread and make a file with the longest isoforms

genes = SeqIO.to_dict(SeqIO.parse(sys.argv[1],"fasta"))
longest = {}

for rec in genes:
	gene = genes[rec].description.split("=")[-1]
	if (gene not in longest) or (longest[gene][0] < len(rec)):
		longest[gene] = (len(genes[rec]),genes[rec].id)

for gene in longest:
	description = genes[longest[gene][1]].id
	genes[longest[gene][1]].id =  gene
	genes[longest[gene][1]].description = description
	print(genes[longest[gene][1]].format("fasta"), end="")
