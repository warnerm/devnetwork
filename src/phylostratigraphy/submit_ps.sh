#!/bin/bash
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 20-00:00
#SBATCH -n 1
snakemake -s Snakefile