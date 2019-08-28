#!/bin/bash
if [ $# -ne 2 ]
then
        echo "Usage: ./DISOPRED <input_file> <out_dir> "
        exit
fi

# $1 the sequence name with suffix .seq
# $2 out directory

RaptorX_HOME=/home/spaenigs/PycharmProjects/encoding_benchmark/snakemake/apps/RaptorX
###### BY TINA, May 5, 2003
BLASTPGP=$RaptorX_HOME/util/BLAST/bin/blastpgp
MAKEMAT=$RaptorX_HOME/util/BLAST/bin/makemat
BLASTDB=$RaptorX_HOME/databases/NR_new/nr90
PSIPREDDIR=$RaptorX_HOME/util/DISOPRED

###### BY TINA, May 5, 2003

DESTDIR=$2

###### BY TINA, May 6, 2003
# first make PSP directory if it doesn't exist
if [ ! -d $DESTDIR ] ; then
    mkdir $DESTDIR
fi
###### BY TINA, May 6, 2003

fulnam=`basename $1`
bname=${fulnam%.*}
rootname=R$bname

echo $rootname 

cat  $1 > $rootname.seq

echo "Running PSI-BLAST with sequence" $1 "...."
#$BLASTPGP -b 0 -j 3 -h 0.0008 -d $BLASTDB -i $1 -C $rootname.chk
#####$BLASTPGP -F T -b 0 -a 22 -j 3 -h 0.001 -d $BLASTDB -i $1 -C $rootname.chk

cp $RaptorX_HOME/tmp/$bname"_sse".chk $rootname.chk
echo $rootname.chk > $rootname.pn
echo $1 > $rootname.sn


#echo $rootname.chk > $rootname.pn
#echo $1 > $rootname.sn
$MAKEMAT -P $rootname

echo Pass1 ...
echo Pass2 ...
$PSIPREDDIR/bin/disopred $rootname $rootname.mtx $PSIPREDDIR/data/


echo "Final output files:" $rootname.diso $rootname.horiz_d
mv $rootname.diso $DESTDIR/$bname.diso
mv $rootname.horiz_d $DESTDIR/$bname.horiz_d

#remove temporary files
echo Cleaning up ....
#-rm -f $rootname.pn $rootname.sn $rootname.aux error.log $rootname.mtx
#-rm -f $rootname.chk 
rm -f $rootname.*
