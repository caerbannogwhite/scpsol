#!/bin/bash
# Parameter sweep: --cut-frequency x --balas-frequency
# Usage: bash benchmarks/sweep.sh <instance_glob> <time_limit> <output_csv>

SOLVER="build/Release/scpsol.exe"
INSTANCES="$1"
TLIMIT="$2"
OUTFILE="$3"

FREQS="0 0.2 0.4 0.6 0.8 1.0"

echo "instance,cut_freq,balas_freq,status,primal,dual,gap,nodes,time" > "$OUTFILE"

for f in $INSTANCES; do
    name=$(basename "$f" .txt)
    for cf in $FREQS; do
        for bf in $FREQS; do
            result=$("$SOLVER" "$f" --time-limit "$TLIMIT" --verbosity 0 \
                --cut-frequency "$cf" --balas-frequency "$bf" 2>&1)
            status=$(echo "$result" | sed -n 's/^Status:[[:space:]]*//p' | tr -d '\r')
            primal=$(echo "$result" | sed -n 's/^Primal:[[:space:]]*//p' | tr -d '\r')
            dual=$(echo "$result" | sed -n 's/^Dual:[[:space:]]*//p' | tr -d '\r')
            gap=$(echo "$result" | sed -n 's/^MIP gap:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')
            nodes=$(echo "$result" | sed -n 's/^Nodes:[[:space:]]*//p' | tr -d '\r')
            wtime=$(echo "$result" | sed -n 's/^Wall time:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')
            echo "$name,$cf,$bf,$status,$primal,$dual,${gap:-0},$nodes,$wtime"
            echo "$name,$cf,$bf,$status,$primal,$dual,${gap:-0},$nodes,$wtime" >> "$OUTFILE"
        done
    done
done
