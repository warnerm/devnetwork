#!/usr/bin/python

#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 1

import pandas as pd
import numpy as np
import sys, getopt

def InOut(argv):
    inputfile=''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'IndivNet_buildup.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i","--ifile"):
            inputfile = arg
        elif opt in ("-o","--ofile"):
            outputfile = arg
    return inputfile, outputfile

#Get ranks for each gene
def getRanks(i):
    d = df.iloc[i]
    s = sorted(range(len(d)), key = lambda k: -d[k])
    return s[1:100] #First value is the index

#Returns the top three genes connected to the given gene
def NextThreeGenes(new_cluster,gene):
    close_genes = np.array(Ranks.iloc[gene, :start_d])
    if sum(gene in new_cluster for gene in close_genes) == 3:
        return new_cluster, []
    else:
        newGenes = [gene for gene in close_genes if not gene in new_cluster]
        new_cluster = np.append(new_cluster, newGenes)
        return new_cluster,newGenes

#Greedily incorporates genes into clusters until no more connections can be made
def AddGenes(new_cluster,gene):
    nextGenes = []
    while True:
        new_cluster, newGenes = NextThreeGenes(new_cluster,gene)
        if len(newGenes) == 0:
            if len(nextGenes) == 0:
                return new_cluster
            else:
                while True:
                    gene = nextGenes[0]
                    nextGenes = nextGenes[1:]
                    if gene not in new_cluster:
                        break
                    elif len(nextGenes) == 0:
                        return new_cluster
        else:
            gene = newGenes[0]
            next = newGenes[1:len(newGenes)]
            nextGenes = np.append(nextGenes,np.array(next))
            nextGenes = nextGenes.astype(int)


#Take the first available gene and start making clusters based on it's top three connections, etc
def makeSmallCluster(genes_remaining):
    start_gene = genes_remaining[0]
    close_genes = np.array(Ranks.iloc[start_gene,:start_d])
    new_cluster = np.append(start_gene,close_genes)
    for gene in close_genes:
        #Add genes until the cluster can't get larger
        new_cluster = AddGenes(new_cluster,gene)

    #Add to existing cluster if applicable
    if sum(gene in genes_remaining for gene in new_cluster) < len(new_cluster):
        for i in range(1,clusterID):
            genes_to_add = np.array(df.Gene[(df['Cluster']==i) & (df['Cluster']!=0)])

            if sum([gene1 in new_cluster for gene1 in genes_to_add]) > 0:
                df.loc[genes_to_add,'Cluster']=clusterID
    #Update to new cluster ID
    df.loc[new_cluster,'Cluster']=clusterID
    genes_remaining = [gene for gene in genes_remaining if not gene in new_cluster]
    return genes_remaining


def ConstructInitialClusters():
    global clusterID
    genes_remaining = np.array(range(Ranks.shape[0]))
    while len(genes_remaining) > 0:
        genes_remaining = makeSmallCluster(genes_remaining)
        clusterID = clusterID + 1

def CombineClusters():
    global clusterID

    #Eventually will update all clusters
    clusters_remaining = np.unique(df.Cluster[df['Cluster']!=0])
    while len(clusters_remaining) > 0:
        current = clusters_remaining[0]

        #Remove cluster form list
        clusters_remaining = np.delete(clusters_remaining,0)

        #Find genes in the cluster
        current_genes = np.array(df.Gene[df['Cluster']==current])
        # print current_genes
        # print d
        # print Ranks.iloc[current_genes,(d-1)]
        #Find the genes reached by incrementing d
        next_genes = np.unique(np.array(Ranks.iloc[current_genes,(d-1)]))
        new_cluster = np.unique(np.append(current_genes,next_genes))

        #Check if new genes are added
        if (len(new_cluster) > len(current_genes)):
            # scan remaining clusters for where those genes are
            for cluster in clusters_remaining:

                # Returns yes if they overlap
                if sum([gene in df.Gene[df['Cluster'] == cluster] for gene in next_genes]) > 0:
                    new_cluster = np.append(new_cluster, df.Gene[df['Cluster'] == cluster])

                    # Remove cluster from list
                    clusters_remaining = [clust for clust in clusters_remaining if clust != cluster]

            # Edit gene cluster ID
            df.loc[new_cluster, 'Cluster'] = clusterID
            clusters_remaining = np.append(clusterID,clusterID)
            clusterID = clusterID + 1

    #Check if giant graph
    if len(np.unique(df['Cluster'])) == 1:
        return True

    else:
        return False

def main(argv):
    inputfile, outputfile = InOut(argv)
    # inputfile = "~/Data/devnetwork/antsTESTpCor.csv"
    # outputfile = "~/Data/devnetwork/antTest.txt"

    global start_d
    start_d = 3 #Start building small graphs with connections of three genes at a time

    #Read in dataframe of pearson correlations
    global df, Ranks, d
    global clusterID
    clusterID = 1 #Iterate as we make clusters
    df = pd.read_csv(inputfile)
    df = abs(df) #Make networks based on magnitude of correlation ('unsigned')
    Ranks = pd.DataFrame([getRanks(i) for i in range(df.shape[0])])
    df = pd.DataFrame(data={'Gene':range(Ranks.shape[0]),'Cluster':np.repeat(0,Ranks.shape[0])})
    ConstructInitialClusters()

    Giant = False
    d = start_d
    while not Giant:
        d = d+1

        #Add a new connection to each gene and check if it forms a giant graph
        Giant = CombineClusters()

    #From calculated d, return the network
    results = {}
    for gene in range(Ranks.shape[0]):
        results[gene] = pd.DataFrame(data = {'Gene1':np.repeat(gene,d),'Gene2':Ranks.iloc[gene,:d]})

    network = pd.concat(results)

    print "final d = " + str(d)
    network.to_csv(outputfile,sep="\t",index=None,header=False)

if __name__ == "__main__":
    main(sys.argv[1:])
    #main("2")