#!/bin/bash
#PBS -l walltime=1:00:00,select=1:ncpus=2:ngpus=1:mem=8gb
#PBS -N batch
#PBS -A st-schieb-1-gpu
 
################################################################################

export NUMBA_CACHE_DIR=$PBS_O_WORKDIR

source ~/.bashrc
module load gcc/9.1.0
conda activate sdecouplings
cd $PBS_O_WORKDIR

SRAND=$RANDOM

python batch.py --adata /home/szhang99/st-schieb-1/sde_couplings/data_repr.h5ad --cellsets /home/szhang99/st-schieb-1/zsteve/reprogramming_batch/cell_sets.gmt --outfile $SRAND.out --srand $SRAND
