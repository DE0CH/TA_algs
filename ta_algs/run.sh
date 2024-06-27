#!/bin/bash -e

FILE=$1
FIRST_LINE=$(head -n 1 $FILE)
D=$(echo $FIRST_LINE | cut -d ' ' -f 1)
N=$(echo $FIRST_LINE | cut -d ' ' -f 2)
GD=$(echo $FIRST_LINE | cut -d ' ' -f 3)
echo "N=$N, D=$D"


make -C ../ 
time ../TA_improved_bardelta -iter 100000 $D $N <(tail -n +2 $FILE)
time ../TA_improved_delta -iter 100000 $D $N <(tail -n +2 $FILE)
