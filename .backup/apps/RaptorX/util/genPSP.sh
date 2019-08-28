#!/bin/bash

# $1 is the sequence name
# $2 is the out directory

if [ $# -lt 2 ]
then
	echo "Usage: ./GenPSP.LINUX <input_file> <out_dir> [cpu_num]"
	exit
fi

RaptorX_CPU=1
if [ $# -gt 2 ]
then
        RaptorX_CPU=$3
fi

RaptorX_HOME=/home/spaenigs/PycharmProjects/encoding_benchmark/snakemake/apps/RaptorX
###### BY TINA, May 5, 2003
BLASTPGP=$RaptorX_HOME/util/BLAST/bin/blastpgp
BLASTDB=$RaptorX_HOME/databases/NR_new/nr90
#BLASTDB=$RAPTOR_HOME/RefSeq/refseq_protein
###### BY TINA, May 5, 2003


###### BY TINA, May 6, 2003
# first make PSP directory if it doesn't exist
if [ ! -d $2 ] ; then
    mkdir $2
fi
###### BY TINA, May 6, 2003

bname=`basename $1 .seq`
$BLASTPGP -a $RaptorX_CPU -b 0 -j 3 -h 0.001 -d $BLASTDB -i $1 -C $2/$bname"_sse".chk 
$BLASTPGP -a $RaptorX_CPU -b 0 -j 3 -h 0.001 -d $BLASTDB -i $1 -R $2/$bname"_sse".chk -C $2/$bname.chk -Q $2/$bname.psp
####$BLASTPGP -a 7 -F T -b 0 -j 5 -h 0.001 -d $BLASTDB -i $1 -C $2/$bname.chk -A 22 -Q $2/$bname.psp
