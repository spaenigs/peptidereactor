#!/usr/bin/env bash

usage()
{
    echo "usage: sysinfo_page [[[-f file ] [-i]] | [-h]]"
}

preprocessing=0
encoding=0

while [ "$1" != "" ]; do
    case $1 in
        -ds | --dataset )         shift
                                  dataset=$1
                                  ;;
        -pt | --part )            shift
                                  part=$1
                                  ;;
        -j | --jobs )             shift
                                  jobs=$1
                                  ;;
        -np | --dry-run )         dryrun=-np
                                  ;;
        -q | --quiet )            quiet=--quiet
                                  ;;
        -pre | --preprocessing )  preprocessing=1
                                  ;;
        -enc | --encoding )       encoding=1
                                  ;;
        -h | --help )             usage
                                  exit
                                  ;;
        * )                       usage
                                  exit 1
    esac
    shift
done

if [ $preprocessing -eq 1 ]; then
  snakemake 00_data/out/$dataset/$dataset\_$part/joblib/$dataset\_$part\_annotated.joblib \
            00_data/out/$dataset/$dataset\_$part/joblib/$dataset\_$part\_annotated_msa.joblib \
	    --config dataset=$dataset part=$part normalize=yes \
            --jobs $jobs \
            --cluster-config cluster.yaml \
            --cluster "qsub -N snmk -S /bin/bash -l h_vmem={cluster.memory} -pe smp {cluster.cores} -l h_rt=172800 -e /scratch/spaenigs/error/ -o /scratch/spaenigs/out" \
            --jobscript scripts/jobscript.sh --latency-wait 60 \
            $dryrun $quiet
fi

if [ $encoding -eq 1 ]; then
  qsub -N scripts/submit.sh
fi

