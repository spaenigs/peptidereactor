#!/usr/bin/env bash

#$ -S /bin/bash
#$ -l h_vmem=500M
#$ -pe smp 32
#$ -l h_rt=259200
#$ -e /scratch/spaenigs/error
#$ -o /scratch/spaenigs/out
#$ -wd /scratch/spaenigs/eb/

. /etc/profile.d/modules.sh
module load tools/R-3.4.3

source /home/spaenigs/.bashrc
conda activate encoding_benchmark

snakemake --config dataset=neuropeptides part=ds1 normalize=yes --cores 32 -np

