#!/usr/local/bin/python2.7
#File takes in list of drosophila developmental proteins and converts them to the OGG names
import pandas as pd
import sys, getopt
import numpy as np


def InOut(argv):
    inputfile=''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:m:o:s:",["ifile=","mfile=","sfile=","ofile="])
    except getopt.GetoptError:
        print 'dmelDev.py -i <inputfile> -m <dmel_map> -s <OGG_map> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i","--ifile"):
            inputfile = arg
        if opt in ("-m", "--mfile"):
            dmel_map = arg
        if opt in ("-s", "--sfile"):
            ogg_map = arg
        elif opt in ("-o","--ofile"):
            outputfile = arg
    return inputfile, dmel_map,ogg_map,outputfile

def main(argv):
    inputfile, dmel_map, ogg_map, outputfile = InOut(argv)
    df = pd.read_csv(inputfile,sep="\t")
    dmel = pd.read_csv(dmel_map,sep="\t",header=None)
    ogg = pd.read_csv(ogg_map,sep=",")
    df.columns = ["Entry","Gene"]
    dmel.columns = ["Entry","Gene"]
    df2 = df.merge(dmel,on="Gene")
    df3 = df2.merge(ogg,left_on="Entry_y",right_on="gene_Dmel_FLYBASE")
    DevOGG = np.unique(np.array(df3['OGG']))
    np.savetxt(outputfile,DevOGG,fmt='%5s',delimiter=',')

if __name__ == "__main__":
    main(sys.argv[1:])


