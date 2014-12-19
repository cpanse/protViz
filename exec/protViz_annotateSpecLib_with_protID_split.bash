#!/bin/bash


# Christian Panse <cp@fgcz.ethz.ch>
# 2014-03-13
# usage:
#   protViz_annotateSpecLib_with_protID_split.bash <speclibfile> <nCPU> <FASTA>


test $# -eq 3 || exit 1
test -x `which protViz_annotateSpecLib_with_protID.py` || { echo "'protViz_annotateSpecLib_with_protID.py' can not be found."; exit 1; }  

SCRATCH=/dev/shm/$$ && mkdir $SCRATCH || exit 1

trap "{ echo 'trapped by return code $? (cleaning up) ...'; rm -vf $SCRATCH/*; rmdir $SCRATCH; exit $?; }" EXIT


function main() {

    input=$1
    nJob=$2
    fasta=$3


    test -s $input || exit 1
    test -s $fasta || exit 1
    test $nJob -lt 1 -o $nJob -gt 50 && nJob=1

    n=`awk -vn=$nJob 'END{print int(NR/n)+1}' $input`

    split $1 -l $n $SCRATCH/speclib. 

    for j in `ls $SCRATCH`;
    do
        ( protViz_annotateSpecLib_with_protID.py --fasta $fasta --speclib $SCRATCH/$j --output $SCRATCH/$j.output ) &
        echo $j
    done

    wait

    sleep 1
    cat $SCRATCH/*.output
}

START_TIME=$SECONDS

main $1 $2 $3 > $1.protein.txt

ELAPSED_TIME=$(($SECONDS - $START_TIME))

echo "$2 $ELAPSED_TIME"

exit 0
