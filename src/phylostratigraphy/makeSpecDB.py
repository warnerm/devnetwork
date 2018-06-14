#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 12

import sqlite3
import re
import sys

spec_names = sys.argv[1]
db = sys.argv[2]
node_file = "nodes.dmp"
name_file = "names_filtered.dmp"


connection = sqlite3.connect(db)

cursor = connection.cursor()

sql_command = """
CREATE TABLE species (
species_number INTEGER PRIMARY KEY,
species VARCHAR(30),
node INT,
parent INT
);"""

cursor.execute(sql_command)

def insert_species(spec):
	spec = spec.strip()
	with open(name_file) as names:
		for line in names:
			if spec in line:
				tax_id = int(re.match('[0-9]+',line).group(0))
				skip = 1
				break;
	if skip == 0: #Didn't find species
		return 0
	with open(node_file) as nodes:
		for line in nodes:
			if int(re.match(r'(\d+)',line).group(0)) == tax_id:
				parent = int(re.search(r"(\d+).*?(\d+)",line).group(2)) #parent node is the second column
				break;
	format_str = """INSERT INTO species (species_number,species,node,parent)
		VALUES (NULL,"{spec}","{tax_id}","{parent}");"""
	sql_command = format_str.format(spec=spec,tax_id=tax_id,parent=parent)
	cursor.execute(sql_command)

with open(spec_names) as input:
	for spec in input:
		insert_species(spec)
	#Parallel(n_jobs = int(threads))(delayed(insert_species) (spec) for spec in input)

connection.commit()

connection.close()

