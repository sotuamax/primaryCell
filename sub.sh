#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=60g
#SBATCH --job-name=nAEC
#SBATCH --time=48:00:00
#SBATCH --mail-user meiyuan.ji@nih.gov 
#SBATCH --mail-type END,FAIL
##SBATCH --gres=lscratch:20

# source myconda 
# mamba activate bio 

ml samtools 

source myconda 
mamba activate bio 

cd /data/jim4/tools 
python3 QC_DHS_temp.py 