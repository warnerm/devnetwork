#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 12

import sys
from Bio import Entrez, SeqIO
import re
import linecache
from joblib import Parallel, delayed
import multiprocessing

name_file = "names_filtered.dmp"
node_file = "nodes.dmp"
cellular_node = 131567 #The node that is cellular organism

def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid"""
    with open(name_file) as input:
        for line in input:
            #print line
            if line.find(species) != -1:
                tax_id = re.match('[0-9]+',line).group(0)
                break;
    return tax_id

def get_tax_data(species):
    print species
    species = species.strip()
    """once we have the taxid, we can fetch the record"""
    current_node = int(get_tax_id(species))
    lineage = ""
    #Iterate through nodes_file until we get to cellular organisms
    while True:
        with open(node_file) as input:
            for line in input:
                #print int(re.match(r'(\d+)',line).group(0))
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
    with open(outfile,'a') as out:
        out.write(species+'\t'+lineage+'\n')
    

spec_names = sys.argv[1]
outfile = sys.argv[2]
threads = 12

with open(spec_names) as input:
    Parallel(n_jobs = int(threads))(delayed(get_tax_data) (spec) for spec in input)






