#!/usr/local/bin/python2.7
import numpy as np, sys, getopt, pandas as pd

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
	        elif opt in ("-o","--ofile"):
	                outfile = arg    
	        elif opt in ("-j","--jfile"):
	               	in2 = arg                    
	return in1,in2,outfile

def main(argv):
	in1,in2,outfile = InOut(argv)
	df1 = pd.read_table(in1,sep="\t")
	df2 = pd.read_table(in2,sep="\t")
	df1.columns=['gene','protein']
	df2.columns=['OGG','protein']
	df3 = pd.merge(df1,df2,how='inner')
	df3.to_csv(outfile,sep=" ",index=None)

if __name__ == "__main__":
	main(sys.argv[1:])