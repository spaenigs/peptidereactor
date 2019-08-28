#!/bin/bash
if [ $# -ne 1 ]
then
        echo "Usage: ./genACC <tgt_name> "
        exit
fi

RaptorX_HOME=/home/spaenigs/PycharmProjects/encoding_benchmark/snakemake/apps/RaptorX
ACCPred=$RaptorX_HOME/util/ACC_Predict/acc_pred
$ACCPred $RaptorX_HOME/tmp/$1.hhm $RaptorX_HOME/tmp/$1.ss2 $RaptorX_HOME/tmp/$1.ss8 $RaptorX_HOME/util/ACC_Predict/model.accpred $RaptorX_HOME/tmp/$1.acc
