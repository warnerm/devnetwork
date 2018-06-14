#Find transcripts that didn't map to anything
import pandas as pd
import sys

genome = sys.argv[1] #list of transcripts
blastMap = sys.argv[2] #results of blasting to ODB
oggMap = sys.argv[3] #ogg - ODBgene map
missing = sys.argv[4]

blast = pd.read_table(blastMap,header=None)
blast.columns = ['transcript','gene']

blast['gene'] = blast['gene'].str.replace('gnl.* ','')
ogg = pd.read_table(oggMap,header=None)
ogg.columns = ['gene','ogg']

together = blast.merge(ogg,on='gene',how='left')
result = together['transcript'][together['ogg'].notnull()].unique()

out = open(missing,'w')

with open(genome) as g:
	for gene in g:
		if gene.strip() not in result:
			out.write(gene)
