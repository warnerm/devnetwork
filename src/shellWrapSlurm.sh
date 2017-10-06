#!/bin/bash
for rep in {1..1000}
do 
	sbatch src/slurmParallel.sh 
	echo rep
done