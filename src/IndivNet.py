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
from scipy import stats as sci

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

#Get ordered list of indices of the genes with the highest sum of pearson correlations
def getTopGenes():
    rows = df.sum()
    s = sorted(range(len(rows)), key=lambda k: -rows[k])
    return s

#Get ranks for each gene
def getRanks(i):
    d = df.iloc[i]
    s = sorted(range(len(d)), key = lambda k: -d[k])
    return s

#Add genes that are connected to the test gene
def newGenes(genes,test,d):
    for g in range(Ranks.shape[0]):

        #skip if already there
        if g in genes: continue

        #Add it if it's connected to the test gene
        elif test in Ranks.iloc[g,:d]: genes =  np.append(genes,g)

    return genes

#Add genes to vector
def addGenes(genes,d,next):


    genes = newGenes(genes,next,d)

    if len(genes) == df.shape[0]:
        bool = True
    else:
        bool = False

    return bool, genes

#Get the next gene to try
def nextGene(genes,topG):
    # Find next gene to add
    ng = 0
    next = topG[0]
    if not next in genes:
        ng = ng + 1
        while (True):
            next = topG[ng]
            if next in genes:
                break
            ng = ng + 1

            #We've run out of options
            if ng == len(topG):
                return next,topG,False

    del topG[ng]
    return next,topG,True

def addNext(genes,next,d):
    genes = np.unique(np.append(genes,Ranks.iloc[next,:d]))
    return genes

#Check if we have a giant network
def CheckIfNetGiant(d):

    #Start with most connected gene each time
    genes = np.array(Ranks.iloc[TopGenes[0],:d])

    #remove used genes from TopGenes
    topG = TopGenes[1:]

    #Add genes to the network, starting with the next highest connected gene that is in the list already
    while(True):
        next,topG,cont = nextGene(genes, topG)
        genes = addNext(genes,next,d)
        #cont returns false if we've exhausted our options (so graph isn't connected)
        if cont:
            bool,  genes = addGenes(genes, d, next)
            if bool: break
        else: return False

    return True

#Print connection list
def getConns(i,d):
    k = Ranks.iloc[i,:d]
    self = np.repeat(i,d)
    return pd.DataFrame({'x':self,'y':k})

#Print network
def getNet(d):
    net = [getConns(i,d) for i in range(Ranks.shape[0])]
    return pd.concat(net)

#Normalize fpkm using hyperbolic sine transformation
def main(argv):
    #inputfile, outputfile = InOut(argv)
    inputfile = "~/Data/devnetwork/beespCor.csv"
    outputfile = "~/Data/devnetwork/beesOne.txt"
    print 's'
    #Start d at 10 just to guess, and go up or down based on whether or not it's a GiantGraph
    d = 10 #variable specifying number of genes each gene is connected to

    #Read in dataframe of pearson correlations
    global df, Ranks, TopGenes
    df = pd.read_csv(inputfile)
    df = abs(df) #Make networks based on magnitude of correlation ('unsigned')
    print 'here'
    Ranks = pd.DataFrame([getRanks(i) for i in range(df.shape[0])])
    print 'here'
    TopGenes = getTopGenes()

    Giant = True
    while Giant:
        print d
        Giant = CheckIfNetGiant(d)
        d = d-1

    #Get the network two steps ago
    net = getNet(d+2)
    print "final" + str(d+2)
    net.to_csv(outputfile,sep="\t",index=None,header=False)

if __name__ == "__main__":
    main("2")
