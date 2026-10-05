#!/bin/bash
# Run all benchmark instances (sets 4, 5, A, B) with default settings.
# Produces a CSV comparable to baseline-no-cuts-no-balas.csv.
#
# Usage: bash benchmarks/run-all.sh [time_limit]
#   time_limit: per-instance time limit in seconds (default: 1800 = 30 min)

SOLVER="build/Release/scpsol.exe"
TLIMIT="${1:-1800}"
OUTFILE="benchmarks/results-dev.csv"

INSTANCES=$(ls data/scp{4,5}{1,2,3,4,5,6,7,8,9,10}.txt data/scp{a,b}{1,2,3,4,5}.txt 2>/dev/null | sort -t/ -k2 -V)

echo "instance,status,primal,dual,mip_gap_pct,nodes,lp_solves,wall_time_s" > "$OUTFILE"

for f in $INSTANCES; do
    name=$(basename "$f" .txt)
    printf "%-10s " "$name"
    # A run whose output is empty or incomplete (see issue #2) is retried
    # instead of silently producing a blank CSV row.
    result=""
    for attempt in 1 2 3; do
        result=$("$SOLVER" "$f" --time-limit "$TLIMIT" --verbosity 0 2>&1)
        rc=$?
        if [ $rc -eq 0 ] && echo "$result" | grep -q '^Status:'; then
            break
        fi
        echo "  (attempt $attempt: exit code $rc, ${#result} bytes of output; retrying)" >&2
        result=""
    done
    if [ -z "$result" ]; then
        echo "FAILED: no result after 3 attempts"
        echo "$name,Failed,,,,,," >> "$OUTFILE"
        continue
    fi
    status=$(echo "$result" | sed -n 's/^Status:[[:space:]]*//p' | tr -d '\r')
    primal=$(echo "$result" | sed -n 's/^Primal:[[:space:]]*//p' | tr -d '\r')
    dual=$(echo "$result" | sed -n 's/^Dual:[[:space:]]*//p' | tr -d '\r')
    gap=$(echo "$result" | sed -n 's/^MIP gap:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')
    nodes=$(echo "$result" | sed -n 's/^Nodes:[[:space:]]*//p' | tr -d '\r')
    lp_solves=$(echo "$result" | sed -n 's/^LP solves:[[:space:]]*//p' | tr -d '\r')
    wtime=$(echo "$result" | sed -n 's/^Wall time:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')
    printf "primal=%-6s dual=%-6s gap=%-10s nodes=%-6s time=%s s\n" \
        "$primal" "$dual" "${gap:-0}%" "$nodes" "$wtime"
    echo "$name,$status,$primal,$dual,${gap:-0},$nodes,$lp_solves,$wtime" >> "$OUTFILE"
done

echo ""
echo "Results written to $OUTFILE"
