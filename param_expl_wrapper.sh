#!/usr/bin/bash

#SBATCH -p medium
#SBATCH --mem=100M
#SBATCH -t 2-0:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -e jmr95_%j.err
#SBATCH -o jmr95_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jromero5@bidmc.harvard.edu

module load gcc/6.2.0
module load conda2/4.2.13
module load python/3.6.0
source activate main_env

inputs[1]=10
inputs[2]=200
inputs[3]=20
inputs[4]=500
inputs[5]=20
inputs[6]=500
inputs[7]=10
inputs[8]=500
inputs[9]=1
inputs[10]=1
inputs[11]=1
inputs[12]=1
inputs[13]=8.5
inputs[14]=0.7
inputs[15]=60
inputs[16]=40
inputs[17]=1000
inputs[18]=3
inputs[19]=-30
inputs[20]=True
inputs[21]=$1
inputs[22]='/n/scratch2/jmr95/parameter_optimization'

srun -c 1 python param_expl.py ${inputs[@]}
