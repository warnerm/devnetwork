import sys
import re, numpy as np
cerana_file = sys.argv[1]
group_file = sys.argv[2]
new_cerana = sys.argv[3]
new_group = sys.argv[4]

originalSeqs = {}

#Make a dictionary of the positions and reference allele
with open(cerana_file,'r') as f:
	for line in f:
		if line[0] == "#":
			continue
		pos = re.split("\t",line)[0] + "_" + re.split("\t",line)[1]
		originalSeqs[pos] = re.split("\t",line)[3]


removeSeqs = [] #Keep track of sequences to remove
removed = 0

out = open(new_group,'w')
out.close()

with open(group_file,'r') as f:
	for line in f:
		if line[0] != "#":
			pos = re.split("\t",line)[0] + "_" + re.split("\t",line)[1]
			if pos in originalSeqs.keys():
				if originalSeqs[pos] != re.split("\t",line)[3]:
					removed = removed+1;
					removeSeqs = np.append(removeSeqs,pos)
					continue
		with open(new_group,'a') as out:
			out.write(line+'\n')

print("Removed "+str(removed)+" sequences from group file")
removed = 0

out = open(new_cerana,'w')
out.close()

with open(cerana_file,'r') as f:
	for line in f:
		if line[0] != "#":
			pos = re.split("\t",line)[0] + "_" + re.split("\t",line)[1]
			if pos in removeSeqs:
				removed = removed+1;
				continue
		with open(new_cerana,'a') as out:
			out.write(line+'\n')

print("Removed "+str(removed)+" sequences from cerana file")
