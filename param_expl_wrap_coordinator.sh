#!/usr/bin/bash

#SBATCH -p short
#SBATCH --mem=10M
#SBATCH -t 0-0:10
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -e jmr95_%j.err
#SBATCH -o jmr95_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jromero5@bidmc.harvard.edu

for ((i=1; i<1000; i=i+1));
do
    sbatch param_expl_wrapper.sh $i.txt
done
