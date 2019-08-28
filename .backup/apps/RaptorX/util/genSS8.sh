#!/bin/bash
if [ $# -ne 1 ]
then
        echo "Usage: ./genSS8 <tgt_name> "
        exit
fi

RaptorX_HOME=/home/spaenigs/PycharmProjects/encoding_benchmark/snakemake/apps/RaptorX
SS8_Pred=$RaptorX_HOME/util/SS8_Predict/bin/run_raptorx-ss8.pl
$SS8_Pred $RaptorX_HOME/tmp/$1.seq -pssm $RaptorX_HOME/tmp/$1.psp -outdir $RaptorX_HOME/tmp/
#mv $1.ss8 $RaptorX_HOME/tmp
#rm $1.*
