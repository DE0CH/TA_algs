#!/bin/bash -e

FILE=$1
FIRST_LINE=$(head -n 1 $FILE)
D=$(echo $FIRST_LINE | cut -d ' ' -f 1)
N=$(echo $FIRST_LINE | cut -d ' ' -f 2)
GD=$(echo $FIRST_LINE | cut -d ' ' -f 3)
echo "N=$N, D=$D"
cargo build --release
export RUST_BACKTRACE=1

make -C ../
cargo run --release --bin=sim_annealing_bardelta -- -i 100000 $D $N <(tail -n +2 $FILE)
# valgrind --leak-check=full --num-callers=20 --tool=memcheck
../TA_improved_bardelta -iter 100000 $D $N <(tail -n +2 $FILE)
../sim_bardelta -iter 100000 $D $N <(tail -n +2 $FILE)
