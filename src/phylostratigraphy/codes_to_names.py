import pandas as pd 
import re
import numpy as np

name_file = "names_filtered.dmp"

def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid"""
    with open(name_file) as input:
        for line in input:
            #print line
            if line.find(species) != -1:
                tax_id = re.match('[0-9]+',line).group(0)
                break;
    return tax_id

infile = "chao_codes.txt"

tbl = pd.read_table(infile,header=None)
taxID = []

for row in range(tbl.shape[0]):
	taxID = np.append(taxID,get_tax_id(tbl.iloc[row,0]))

tbl = pd.concat([tbl,pd.DataFrame(taxID)],axis=1)
tbl.to_csv('chao_codes_edit.txt',sep='\t',index=None,header=None)