#!/bin/bash -e

FILE=$1
FIRST_LINE=$(head -n 1 $FILE)
D=$(echo $FIRST_LINE | cut -d ' ' -f 1)
N=$(echo $FIRST_LINE | cut -d ' ' -f 2)
GD=$(echo $FIRST_LINE | cut -d ' ' -f 3)
echo "N=$N, D=$D"
cargo build --release
export RUST_BACKTRACE=1
time ./target/release/bardelta $D $N <(tail -n +2 $FILE) --seed $RANDOM
time ../TA_improved_bardelta $D $N <(tail -n +2 $FILE)
time ./target/release/delta $D $N <(tail -n +2 $FILE) --seed $RANDOM
time ../TA_improved_delta $D $N <(tail -n +2 $FILE)

echo "ground truth $GD"
