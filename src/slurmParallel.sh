#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --array=1-5
#SBATCH -t 0-24:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mem-per-cpu=4GB
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=warnerm@sas.upenn.edu # send-to addressfor i in {1..100000}; do
Rscript src/analysis_parallelTests.R