from __future__ import division
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import fileinput

#assumes that the gene id is the first line


def changes(gene):
    """count number of silent and replacement sites"""
    degeneracy = {"A" : 4,"R" : 6,"N" : 2,"D" : 2,"C" : 2,"Q" : 2,"E" : 2,"G" : 4,"H" : 2,"I" : 2,"M" : 2,"L" : 6,"K" : 2,"F" : 2,"P" : 4,"S" : 6,"T" : 4,"W" : 1,"Y" : 2,"V" : 4,"*" : 0}
    Trepl = 0
    Tsil = 0
    for pos in range(0,len(gene) - len(gene) % 3 ,3):
        #check for possible transtions
        aa = str(gene.seq[pos:pos+3].translate())
        if aa not in degeneracy:
            continue
        silent  = (degeneracy[aa]-1)/3  #number of silent sites
        Trepl += 3 - silent
        Tsil += silent

    return(Trepl, Tsil)

#read reference data and return Tsil, Trepl for longest isoform
genes = {}
print("gene,isoform,Trepl,Tsil")
for rec in SeqIO.parse(fileinput.input(),"fasta"):
    Trepl, Tsil = changes(rec)
    print('%s,%s,%.4f,%.4f' % (rec.id,rec.description.split(" ")[-1],Trepl, Tsil))
