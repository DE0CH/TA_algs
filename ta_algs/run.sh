#!/bin/bash -e

FILE=$1
FIRST_LINE=$(head -n 1 $FILE)
D=$(echo $FIRST_LINE | cut -d ' ' -f 1)
N=$(echo $FIRST_LINE | cut -d ' ' -f 2)
GD=$(echo $FIRST_LINE | cut -d ' ' -f 3)
echo "N=$N, D=$D"
cargo build --release
make -C ../
cargo flamegraph --release --bin=delta -- $D $N <(tail -n +2 $FILE)
perf record -F 99 -g -- ../TA_improved_delta $D $N <(tail -n +2 $FILE)
perf script > out.perf
../../FlameGraph/stackcollapse-perf.pl out.perf > out.folded
../../FlameGraph/flamegraph.pl out.folded > flamegraph_c.svg
time ./target/release/delta $D $N <(tail -n +2 $FILE)
time ../TA_improved_delta $D $N <(tail -n +2 $FILE)

echo "ground truth $GD"
