#!/bin/bash
# Targeted parameter sweep on scpb* validation set
# Tests only promising (cut_freq, balas_freq) combos from training

SOLVER="build/Release/scpsol.exe"
TLIMIT="$1"
OUTFILE="$2"

echo "instance,cut_freq,balas_freq,status,primal,dual,gap,nodes,time" > "$OUTFILE"

# Top candidates from training (scpa* sorted by time):
# (0, 0), (0, 0.6), (0.8, 1.0), (1.0, 0.4), (0.2, 0.2),
# (0.4, 0), (0, 0.8), (0.2, 0.4), (1.0, 0.2), (0, 1.0)
# Plus current defaults (0.4, 0.2) and some other interesting ones
COMBOS="0,0 0,0.6 0,0.8 0,1.0 0.2,0.2 0.2,0.4 0.4,0 0.4,0.2 0.6,0.4 0.8,0.4 0.8,1.0 1.0,0.2 1.0,0.4"

for f in data/scpb*.txt; do
    name=$(basename "$f" .txt)
    for combo in $COMBOS; do
        cf=$(echo "$combo" | cut -d, -f1)
        bf=$(echo "$combo" | cut -d, -f2)
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
