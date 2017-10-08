#!/bin/bash#
#SBATCH -p compute # partition (queue)
python2.7 ../src/IndivNet.py -i bees.tpm.txt -o bees.net.txt