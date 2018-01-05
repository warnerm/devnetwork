#Generates OGG map in orthoclust framework
import pandas as pd
import numpy as np
ogg = pd.read_table("~/GitHub/devnetwork/data/HymOGG_hym.csv",sep = " ")

#Get 1-1-1 orthologs
ogg = ogg.drop_duplicates('OGG',keep=False)
ogg = ogg.drop_duplicates('gene_Mphar',keep=False)
ogg = ogg.drop_duplicates('gene_Amel',keep=False)

bee = pd.read_table("~/Github/devnetwork/data/bees.tpm.txt")
ant = pd.read_table("~/Github/devnetwork/data/ants.tpm.txt")
bee['bee_ind']=range(0,bee.shape[0])
ant['ant_ind']=range(0,ant.shape[0])

oggB = pd.merge(ogg,bee,how='inner',left_on = 'gene_Amel', right_on = 'id')
oggC = pd.merge(oggB,ant,how='inner',left_on = 'gene_Mphar', right_on = 'id')
oggC['ogg_ind']=range(0,oggC.shape[0])

ret = oggC[['bee_ind','ant_ind']]
ret['sp1'] = np.repeat(1,ret.shape[0])
ret['sp2'] = np.repeat(2,ret.shape[0])
ret['couple'] = np.repeat(1,ret.shape[0])
ret = ret[['sp1','sp2','couple','bee_ind','ant_ind']]

ret.to_csv('~/Data/devnetwork/OGGmap_orthoclust.txt',sep='\t',index=None,header=False)