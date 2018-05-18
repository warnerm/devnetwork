#!/usr/bin/python
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 1

import sys
from Bio import Entrez, SeqIO
import re

#From https://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython
def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(' ', "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    print record['Count']
    if int(record['Count']) == 0:
        return 0
    else:
        return record['IdList'][0]

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)

Entrez.email = ""
output = open('ncbiTest_edit.fa','w')
output.close()

with open('nr_ncbi') as input:
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
                tax = get_tax_data(tax_id)
                lineage = tax[0]['Lineage']
                lineage = re.sub('cellular organisms; ','',lineage)
                line = line + ' | ' + '[' + lineage + ']\n'
                skip = False
            with open('ncbiTest_edit.fa','a') as output:
                output.write(line)
        else:
            if not skip: 
                with open('ncbiTest_edit.fa','a') as output:
                    output.write(line)








