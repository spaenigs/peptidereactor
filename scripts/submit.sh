#!/usr/bin/env bash

#$ -S /bin/bash
#$ -l h_vmem=750M
#$ -pe smp 16
#$ -l h_rt=86400
#$ -e /scratch/spaenigs/error
#$ -o /scratch/spaenigs/out

. /etc/profile.d/modules.sh
module load tools/R-3.4.3

source .bashrc
conda activate encoding_benchmark

snakemake --config dataset=neuropeptides part=ds1 normalize=yes --cores 16 -np

