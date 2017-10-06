#!/bin/bash
for rep in {1..1000}:
do 
	sbatch slurmParallel.sh 
done