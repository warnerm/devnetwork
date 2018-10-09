#takes two codon codon files
#then runs paml in pairwise mode to estimate dn/ds
from Bio.Phylo.PAML import codeml
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from itertools import izip
import os,sys,pdb

nuc1 = sys.argv[1]
nuc2 = sys.argv[2]
with open("scratch/tree.tre","w") as outfile:
    outfile.write("(cerana_,mellifera);\n")
    #assume chinense is first
for rec1,rec2 in izip(SeqIO.parse(sys.argv[1],"fasta"), SeqIO.parse(sys.argv[2],"fasta")):
    moduloRec =len(rec1) % 3
    if moduloRec != 0:
        rec1 = rec1[:-1 * moduloRec ]
        rec2 = rec2[:-1 * moduloRec ]
    gene = rec1.description.split("=")[1]
    if str(rec2[:-3].seq.translate()).find("*") > -1:
        print rec2[:-3].seq.translate()
        #skip genes with in-frame stop codons in chinense, as potentially problematic
        print "{}\tSTOP".format(gene)
        continue
    rec1.id = "cerana_"
    rec2.id = "mellifera"
    AlignIO.write(MultipleSeqAlignment([rec1[:-3],rec2[:-3]]),"scratch/alignment.phy","phylip-sequential")
    os.system("""perl -p -i -e 's/ /\n/' scratch/alignment.phy""") #replace spaces by newlines for codeml
cml = codeml.Codeml()
cml.alignment = "scratch/alignment.phy"
cml.working_dir = "./scratch"
cml.tree = "scratch/tree.tre"
cml.out_file = "scratch/out.txt"
cml.set_options(seqtype=1,
        verbose=1,
        noisy=3,
        model=1,
        runmode=-2,
        Mgene=0,
        NSsites=[0],
        CodonFreq=2,
        cleandata=1)
cml.run(verbose=True)

    # parse codemls
    with open(cml.out_file) as results: 
        for line in results:
            if line.find("dN/dS=") > -1:
                line = line.split()
                print "%s\t%s" % (gene, line[line.index("dN/dS=")+1])
