#!/usr/bin/python

#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 40

#File takes in tpm values and constructs a network based on pearson correlations among genes
#Finds a value for d, the threshold of pearson correlations, so that all genes are in a giant graph
from compiler.ast import nodes

import pandas as pd
import sys, getopt
import numpy as np
from joblib import Parallel, delayed
import multiprocessing

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

#Normalize fpkm using hyperbolic sine transformation
def hypsine ( d ):
    return np.log(d + np.sqrt(d ** 2 + 1))

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

def callConnections(dfCor,d,gene):
    GeneCor = dfCor.iloc[gene]
    order = GeneCor.sort_values(ascending=False)
    keep = order[1:d + 1].index[0:d].tolist()
    g = np.repeat(gene,d)
    ret = pd.DataFrame({'gene':g,'conn':keep})
    return ret

#make network based on absolute value of pearson correlation
def MakeNetwork(dfCor,d):
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(callConnections)(dfCor,d,gene) for gene in range(dfCor.shape[0]))
    results = pd.concat(results)
    return results

def main(argv):
    inputfile, outputfile = InOut(argv)

    #Start d at 10 just to guess, and go up or down based on whether or not it's a GiantGraph
    d = 10 #variable specifying number of genes each gene is connected to

    #Read in dataframe of pearson correlations
    df = pd.read_csv(inputfile)
    df = abs(df) #Make networks based on magnitude of correlation ('unsigned')
    net = MakeNetwork(df, d)
    if (GiantGraph(net)): direction = 0
    else: direction = 1
    while True:
        print d
        net = MakeNetwork(df,d)

        #direction == 1 means we are making bigger graphs
        if (direction):
            if GiantGraph(net):
                break
            else: d = d + 1

        #direction == 0 means we are making smaller graphs
        else:
            if not GiantGraph(net):
                net = MakeNetwork(df,d+1) #We've gone one step to far
            else: d = d - 1

    print "final" + str(d)
    net.to_csv(outputfile,sep="\t",index=None,header=False)

if __name__ == "__main__":
    main(sys.argv[1:])
