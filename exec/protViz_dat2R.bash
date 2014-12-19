#!/bin/bash

INPUTFILE=$1 
MODFILE=/usr/local/mascot/mascot_2_4/config/mod_file

test -f $INPUTFILE || { echo "$INPUTFILE does not exists." ; exit 1; }
test -f $MODFILE || { echo "can not find $MODFILE," ; exit 1; }

MYFILENAME=`basename $INPUTFILE .dat`

test -f `which protViz_mascotDat2RData.pl` \
&& DAT2R=protViz_mascotDat2RData.pl


$DAT2R -d=$INPUTFILE -m=$MODFILE || { echo "$0 failed." ; exit 1; }
