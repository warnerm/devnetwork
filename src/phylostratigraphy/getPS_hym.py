import sys,re
from joblib import Parallel, delayed
import numpy as np
import sqlite3
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
from subprocess import call


global outfile
global main_tax


infile = sys.argv[1]
outfile = sys.argv[2]
threads = sys.argv[3]
blastDB = sys.argv[4]
seqfile = sys.argv[6] #species protein file
global max_length
global l1

main_tax = int(sys.argv[5]) #Species tax_id

#Return tax_id of all blast hits in database
def blastTranscript(transcript,blastDB,seqfile):
	transcript = transcript.strip()
	spec = seqfile.replace('_prot.fa','')
	blastDB = blastDB.replace('.phr','')

	#Make fasta file of the individual protein
	seqiter = SeqIO.parse(open(seqfile),'fasta')
	SeqIO.write((seq for seq in seqiter if seq.id in transcript), "temp"+transcript+".fa", "fasta")

	blastp_cline = NcbiblastpCommandline(query="temp"+transcript+".fa",db=blastDB,evalue=1e-10,
											outfmt=5,out="blast"+transcript+".xml")

	stdout, stderr = blastp_cline()

	result_handle = open("blast"+transcript+".xml")
	blast_record = NCBIXML.read(result_handle)
	E_VALUE_THRESH = 1e-10

	alignments = []

	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESH:
				with open(outfile,'a') as f:
					alignments = np.append(alignments,str(alignment.title))

	call(["rm","temp"+transcript+".fa","blast"+transcript+".xml"])
	alignments2 = [re.split(' ',aln)[1] for aln in alignments]
	alignments3 = [re.split('_',aln)[0] for aln in alignments2]
	return alignments3


#Compare lineages and find the first point of divergence
def compLineage(lineage2):
	ps = max_length
	l2 = re.split('; ',lineage2)
	for i in range(max_length):
		s1 = l1[i]
		s2 = l2[i]
		if s1 != s2:
			return i
	return max_length


#Find the last node at which two groups match based on lineage from SQL db
def getMatch(tax):
	if main_tax == tax or int(tax) == 9478:
		return 40 ##Just return some high number so it doesn't get picked
	else:
		connection = sqlite3.connect('../data/taxonomy.db')
		cursor = connection.cursor()
		format_str = """SELECT lineage FROM taxonomy WHERE tax_id = {tax_id}"""
		sql_command = format_str.format(tax_id=int(tax))
		cursor.execute(sql_command)
		try:
			lineage2 = cursor.fetchone()[0]
		except:
			connection.close()
			return 40
		connection.close()
		return compLineage(lineage2)

#Get phylostrata of most distant match for a given OGG
def getPS(transcript,blastDB,seqfile):
	OGGset = blastTranscript(transcript,blastDB,seqfile)
	OGGset = [codes[name] for name in OGGset] #get tax id
	Matches = [getMatch(tax) for tax in OGGset]
	if len(Matches) < 1:
		ps = 'not_found'
	else:
		ps = np.min(Matches)
	with open(outfile,'a') as out:
		out.write(transcript.strip()+'\t'+str(ps)+'\n')

out = open(outfile,'w')
out.close()


if main_tax == 7227: #Skip the whole process with Drosophila
	with open(outfile,'a') as out:
		out.write('No results for Drosophila melanogaster')

else:
	connection = sqlite3.connect('../data/taxonomy_hym.db')
	cursor = connection.cursor()
	format_str = """SELECT lineage FROM taxonomy_hym WHERE tax_id = {tax_id}"""
	sql_command = format_str.format(tax_id=int(main_tax))
	cursor.execute(sql_command)
	lineage = cursor.fetchone()[0]
	connection.close()

	max_length = len(re.split('; ',lineage))
	l1 = re.split('; ',lineage)

	codes = {}
	with open("../data/chao_codes_edit.txt") as ch:
		for line in ch:
			name = re.split('\t',line)[1]
			tax_id = re.split('\t',line)[2].strip()
			codes[name] = tax_id

	with open(infile) as geneList:
		for gene in geneList:
			Parallel(n_jobs=int(threads))(delayed(getPS)(gene,blastDB,seqfile) for gene in geneList)








