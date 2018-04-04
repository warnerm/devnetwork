from glob import glob
from collections import defaultdict
import sys,pdb
""" 
collects fpkm for either genes or isoforms produced by an RSEM analysis, 
and tabulate expected counts across samples.
INPUT: takes `genes` or `isoforms` as input, and a location where they can be found
OUTPUT: print csv table to stdout
"""

if sys.argv[1] != 'genes' and sys.argv[1] != 'isoforms' or not sys.argv[2]:
	raise ValueError

counts = defaultdict(list)
samples = []
 
for file in sorted(glob(sys.argv[2]+"*"+sys.argv[1]+".results")):
	infile = open(file)
	samples.append(file.split("/")[-1].rstrip("."+sys.argv[1]+".results")) # save sample id
	infile.next()
	for line in infile:
		line = line.split()
		counts[line[0]].append(line[-1])
	infile.close()

print ",".join([sys.argv[1][:-1]+"_id"] + samples)
for transcript in counts:
	print ",".join([transcript] + counts[transcript])



