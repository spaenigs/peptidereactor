#!/usr/bin/env bash

TOKEN=$1
INPUT=$2
OUTPUT=$3
TYPE=$4

RND=$(shuf -i 100000-999999 -n 1)

OLD_WD=$PWD
PATH_TO_NODE=nodes/encodings/electrostatic_hull/;
NEW_WD=data/temp/$TOKEN/$RND/; mkdir $NEW_WD;

cp $INPUT $NEW_WD/template.pqr;
cp $PATH_TO_NODE/config/apbs_$TYPE.input $NEW_WD;
cd $NEW_WD;

apbs apbs_$TYPE.input 1> /dev/null 2> /dev/null;

cp $TYPE.dx $OLD_WD/$OUTPUT;
cd - 1> /dev/null;

rm -r $NEW_WD;