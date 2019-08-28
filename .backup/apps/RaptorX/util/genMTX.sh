#!/bin/bash
if [ $# -ne 1 ]
then
        echo "Usage: ./genMTX <tgt_name> "
        exit
fi

RaptorX_HOME=/home/spaenigs/PycharmProjects/encoding_benchmark/snakemake/apps/RaptorX
BLAST_HOME=$RaptorX_HOME/util/BLAST
MakeMat=$RaptorX_HOME/util/BLAST/bin/makemat

#echo $RaptorX_HOME/tmp/$1.seq > $RaptorX_HOME/tmp/$1.sn
#echo $RaptorX_HOME/tmp/$1.chk > $RaptorX_HOME/tmp/$1.pn
#$MakeMat -P $RaptorX_HOME/tmp/$1
#rm $RaptorX_HOME/tmp/$1.sn
#rm $RaptorX_HOME/tmp/$1.pn


cd $RaptorX_HOME/tmp
echo $1.seq > $1.sn
echo $1.chk > $1.pn
$BLAST_HOME/bin/makemat -P $1
rm $1.sn
rm $1.pn
