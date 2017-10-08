#!/usr/local/bin/python2.7
#File takes in tpm values and constructs a network based on pearson correlations among genes
from compiler.ast import nodes

import pandas as pd
import sys, getopt
import numpy as np


def InOut(argv):
    inputfile=''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'IndivNet.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i","--ifile"):
            inputfile = arg
        elif opt in ("-o","--ofile"):
            outputfile = arg
    return inputfile, outputfile

#Take in dataframe of nodes and check if they are all connected
def GiantGraph(net):
    nodes = []
    for node in net[net.columns[0]].unique(): #iterate through each gene
        length = len(nodes)
        sub = net[net[net.columns[0]] == node]
        for row in range(sub.shape[0]): #Add new nodes to the list
            if sub.iloc[row,1] in nodes:
                continue
            else:
                nodes.append(sub.iloc[row,1])

        newlen = len(nodes)
        if (newlen == length): #If the list hasn't grown, the graph is disconnected
            return False
        else:
            if newlen == len(net[net.columns[0]].unique()): #if the list is the full length, graph is connected
                return True
            else:
                continue



#make network based on absolute value of pearson correlation
def MakeNetwork(dfCor,d):
    net= pd.DataFrame([[0,0]])
    for gene in range(dfCor.shape[0]):
        GeneCor = dfCor.iloc[gene]
        order = GeneCor.sort_values(ascending=False)
        keep = order[1:d+1].index[0:d].tolist() #skip first value because that is self-correlation
        dfSm = pd.DataFrame([[order.index[0]]*d,keep])
        dfSm = dfSm.transpose()
        net=pd.concat([net,dfSm])
    return net[1:]

def main(argv):
    inputfile, outputfile = InOut(argv)
    df = pd.read_table(inputfile, sep="\t")
    filter = df[df.columns[0]].str.contains("ERCC")
    df = df[~filter]
    genes = df[df.columns[0]]
    df = df.drop(df.columns[[0]],axis=1)
    sums = df.sum(axis=1)
    df = df[sums > 0]
    dfCor = np.corrcoef(np.array(df))
    dfCor = pd.DataFrame(dfCor)
    d = 1 #variable specifying number of genes each gene is connected to
    while True:
        net = MakeNetwork(dfCor.abs(),d)
        if (GiantGraph(net)):
            break
        else:
            d = d+1

    net.to_table(outputfile,sep="\t",index=None)

if __name__ == "__main__":
    main(sys.argv[1:])
