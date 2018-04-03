#!/bin/bash
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 10-00:00
#SBATCH -n 1
SRAKEY=`cat samples_jasper.txt | awk '{{print $4}}' | sed '1d' | tr -d ' '`
echo $SRAKEY
for sra_key in $SRAKEY:
do
	fastq-dump $sra_key
done