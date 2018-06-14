import sqlite3
import re
import sys
import os.path

spec_nodes = sys.argv[1]
node_file = "../data/nodes.dmp"
cellular_node = 131567 #The node that is cellular organism
all_node = 1 #Viruses will hit this node

connection = sqlite3.connect('../data/taxonomy_hym.db')

cursor = connection.cursor()

#Create the sql database if it doesn't exist
#if not os.path.isfile("taxonomy.db"):
sql_command = """
CREATE TABLE taxonomy_hym (
species_number INTEGER PRIMARY KEY,
tax_id INT,
lineage VARCHAR(200)
);"""
cursor.execute(sql_command)

connection.commit()

connection.close()

def addTax(node):
	tax_id = node.strip()
	current_node = int(tax_id)
	lineage = ""
	#Iterate through nodes_file until we get to cellular organisms
	while True:
		if current_node == cellular_node or current_node == 1: break;
		found = 0
		with open(node_file) as input:
			for line in input:
				if int(re.match(r'(\d+)',line).group(0)) == current_node:
					current_node = int(re.search(r"(\d+).*?(\d+)",line).group(2)) #parent node is the second column
					lineage = str(current_node) + "; " + lineage
					found = 1
					break;
			if found == 0: return 0
	
	lineage = lineage[:-2] #Remove last semi-colon
	connection = sqlite3.connect('../data/taxonomy_hym.db')

	cursor = connection.cursor()

	format_str = """INSERT INTO taxonomy_hym (species_number,tax_id,lineage)
		VALUES (NULL,"{tax_id}","{lineage}");"""
	sql_command = format_str.format(tax_id=tax_id,lineage=lineage)
	cursor.execute(sql_command)

	connection.commit()

	connection.close()

#Make code to tax_id dictionary
codes = {}
with open("../data/chao_codes_edit.txt") as ch:
	for line in ch:
		name = re.split('\t',line)[1]
		tax_id = re.split('\t',line)[2].strip()
		codes[name] = tax_id

with open(spec_nodes) as spec:
	for line in spec:
		species = line.replace('_edit.fa','')
		addTax(codes[species.strip()])







