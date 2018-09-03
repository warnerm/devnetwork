from Bio import SeqIO
import sys, pdb, sets

fasta=sys.argv[1] # file with longest isoforms
annot=sys.argv[2] # file with snpEff output

#read longest isoform
isoforms = {}
for rec in SeqIO.parse(fasta,"fasta"):
    isoforms[rec.description.split(" ")[1]] = rec.id

def isoform_edit(isoform):
    try:
        isoform = isoform.split("|")
        isoName = isoform[6]
        # pdb.set_trace()
        #ignore non-longest isoforms
        if isoName not in isoforms:
            # pdb.set_trace()
            return
        # use effect to estim
        if isoform[2] == "LOW" or isoform[2] == "NONE":
            effect = "S"
        else:
            effect  = "N"
        print("{},{},{},{},{}".format(chrom,pos,effect,isoforms[isoName],isoName))
    except:
        print "Error"


with open(annot) as annotfile:
    for line in annotfile:
        if line[0] == "#" or line.find("ANN=") == -1:
            continue
        line = line.split("\t")
        chrom, pos, eff  = line[0], line[1], line[7]
        for isoform in eff.split(","):
            isoform_edit(isoform)