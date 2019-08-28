#!/bin/bash
if [ $# -lt 1 ]
then
        echo "Usage: ./genHMM <tgt_name> [cpu_num]"
        exit
fi

RaptorX_CPU=1
if [ $# -gt 1 ]
then
	RaptorX_CPU=$2
fi

RaptorX_HOME=/home/spaenigs/PycharmProjects/encoding_benchmark/snakemake/apps/RaptorX
cd $RaptorX_HOME/util/HHpred
f=`echo $(basename $1 .seq)`
./buildali2.pl $RaptorX_HOME/tmp/$f.seq -cpu $RaptorX_CPU
./hhmake -i $RaptorX_HOME/tmp/$f.a3m
