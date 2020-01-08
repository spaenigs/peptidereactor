#!/usr/bin/env bash

export TOKEN=$1
export INPUT=$2
export OUTPUT=$3
export TYPE=$4

export OLD_WD=$PWD
export PATH_TO_NODE=nodes/encodings/electrostatic_hull/;
export NEW_WD=data/temp/$TOKEN/$RANDOM/; mkdir $NEW_WD;

cp $INPUT $NEW_WD/template.pqr;
cp $PATH_TO_NODE/config/apbs_$TYPE.input $NEW_WD;
cd $NEW_WD;

apbs apbs_$TYPE.input 1> /dev/null 2> /dev/null;

cp $TYPE.dx $OLD_WD/$OUTPUT;
cd - 1> /dev/null;

rm -r $NEW_WD;