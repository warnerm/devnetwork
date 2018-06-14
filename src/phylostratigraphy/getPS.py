import sys,re
from joblib import Parallel, delayed
import numpy as np
import sqlite3

global outfile

infile = sys.argv[1]
outfile = sys.argv[2]
threads = sys.argv[3]
odbDB = sys.argv[4]
global main_tax
global max_length
global l1

# #Get list of all species with the OGG
def getRepresentation(OGG):
	speciesList = []
	with open(odbDB) as odb:
		for line in odb:
			if OGG in line:
				species = re.split('\t',line)[1]
				species = re.split(':',species)[0]
				speciesList = np.append(speciesList,species)
	return speciesList

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
def getPS(OGG):
	OGG = re.split('\t',OGG)[1].strip()
	OGGset = getRepresentation(OGG) #Return list of species with the OGG
	Matches = [getMatch(tax) for tax in OGGset]
	ps = np.min(Matches)
	with open(outfile,'a') as out:
		out.write(OGG+'\t'+str(ps)+'\n')

with open(infile) as geneList:
	first_line = geneList.readline()
	main_tax = re.split(':',first_line)[0]

out = open(outfile,'w')
out.close()

connection = sqlite3.connect('../data/taxonomy.db')
cursor = connection.cursor()
format_str = """SELECT lineage FROM taxonomy WHERE tax_id = {tax_id}"""
sql_command = format_str.format(tax_id=int(main_tax))
cursor.execute(sql_command)
lineage = cursor.fetchone()[0]
connection.close()

max_length = len(re.split('; ',lineage))
l1 = re.split('; ',lineage)

with open(infile) as geneList:
	Parallel(n_jobs=int(threads))(delayed(getPS)(gene) for gene in geneList)
