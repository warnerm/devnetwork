#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 1

import sys
from Bio import Entrez, SeqIO
import re
import linecache

name_file = "names_filtered.dmp"
node_file = "nodes.dmp"
cellular_node = 131567 #The node that is cellular organism

def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid"""
    with open(name_file) as input:
        for line in input:
            if species in line:
                tax_id = re.match('[0-9]+',line).group(0)
    return tax_id

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    current_node = int(taxid)
    lineage = ""
    #Iterate through nodes_file until we get to cellular organisms
    while True:
        print current_node
        with open(node_file) as input:
            for line in input:
                #print line
                #print re.match(r'(\d+)',line).group(0)
                if int(re.match(r'(\d+)',line).group(0)) == current_node:
                    current_node = int(re.search(r"(\d+).*?(\d+)",line).group(2)) #parent node is the second column
                    break;
        if current_node == cellular_node: break;
        with open(name_file) as input:
            for line in input:
                if int(re.match(r'(\d+)',line).group(0)) == current_node:
                    lineage = re.search('[A-Za-z]+',line).group(0) + "; " + lineage
                    break; #Necessary because file contains synonyms after the first level
    lineage = lineage[:-2] #Remove last semi-colon
    return lineage

with open("ncbiTest.fa") as input:
    for line in input:
        if line.startswith('>'):
            line = re.sub(' .*\[',' | [',line).strip()
            species = re.sub('.*\[','',line)
            species = re.sub('\]','',species)
            tax_id = get_tax_id(species)
            if tax_id == 0:
                line = ''
                skip = True
            else:
                lineage = get_tax_data(tax_id)
                line = line + ' | ' + '[' + lineage + ']\n'
                skip = False
            with open('ncbiTest_edit.fa','a') as output:
                output.write(line)
        else:
            if not skip: 
                with open('ncbiTest_edit.fa','a') as output:
                    output.write(line)            
    






