#!/usr/local/bin/python2.7
import numpy as np, sys, getopt, pandas as pd
import re

def InOut(argv):
	in1=''
	in2=''
	outfile = ''
	try:
	        opts, args = getopt.getopt(argv,"h:i:j:o:",["ifile=","jfile=","ofile="])
	except getopt.GetoptError:
	        print 'combineOGGmaps.py -i <input1> -j <input2> -o <outfile>'
	        sys.exit(2)
	for opt, arg in opts:
	        if opt in ("-i","--ifile"):
	                in1 = arg 
	        elif opt in ("-j","--jfile"):
	                in2 = arg                  
	        elif opt in ("-o","--ofile"):
	                outfile = arg            
	return in1,in2,outfile

def main(argv):
	in1,in2,outfile = InOut(argv)
	df1 = pd.read_csv(in1,sep=" ")
	df2 = pd.read_csv(in2,sep=" ")
	df1 = df1.drop('protein',axis=1)
	df2 = df2.drop('protein',axis=1)
	name1=re.search('OGG_(\w+)_',in1)
	name2=re.search('OGG_(\w+)_',in2)
	df1.columns=['gene_'+name1.group(1),'OGG']
	df2.columns=['gene_'+name2.group(1),'OGG']
	df3 = pd.merge(df1,df2,how='inner')
	df3 = df3.drop_duplicates()
	cols = df3.columns.tolist()
	cols = cols[1:] + cols[:1]
	df3 = df3[cols]
	df3.to_csv(outfile,sep=" ",index=None)

if __name__ == "__main__":
	main(sys.argv[1:])