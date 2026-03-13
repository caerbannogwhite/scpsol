#!/bin/bash
# Overnight parameter sweep with parallelism.
# Training: scpb1-5, scpnre1-3, scpnrf1-3 (11 instances)
# Validation: scpnre4-5, scpnrf4-5 (4 instances)
#
# Usage: bash benchmarks/sweep-overnight.sh [train|validate] [time_limit] [parallel]
#   train:     run training sweep (default)
#   validate:  run validation on specified combos
#   time_limit: per-instance time limit in seconds (default 600)
#   parallel:  max parallel jobs (default 6)

SOLVER="build/Release/scpsol.exe"
MODE="${1:-train}"
TLIMIT="${2:-300}"
PARALLEL="${3:-6}"

if [ "$MODE" = "train" ]; then
    INSTANCES="data/scpb1.txt data/scpb2.txt data/scpb3.txt data/scpb4.txt data/scpb5.txt data/scpnre1.txt data/scpnre2.txt data/scpnre3.txt data/scpnrf1.txt data/scpnrf2.txt data/scpnrf3.txt"
    OUTFILE="benchmarks/sweep-overnight-train.csv"
elif [ "$MODE" = "validate" ]; then
    INSTANCES="data/scpnre4.txt data/scpnre5.txt data/scpnrf4.txt data/scpnrf5.txt"
    OUTFILE="benchmarks/sweep-overnight-validate.csv"
else
    echo "Usage: $0 [train|validate] [time_limit] [parallel]"
    exit 1
fi

# Build the list of (instance, params) jobs
JOBFILE=$(mktemp)

# Parameter grid:
#   cut_freq: 0, 0.2
#   balas_freq: 0, 0.6
#   branch: most_fractional (eta/sb ignored), reliability (eta x sb grid)
#   eta: 2, 4, 8
#   sb: 4, 8, 16
# Total combos per instance: 2 * 2 * (1 + 3*3) = 40
# Total runs (train): 40 * 11 = 440
# At 600s/run, 6 parallel: ~12h

for f in $INSTANCES; do
    for cf in 0 0.2; do
        for bf in 0 0.6; do
            echo "$f most_fractional $cf $bf 0 0" >> "$JOBFILE"
            for eta in 2 4 8; do
                for sb in 4 8 16; do
                    echo "$f reliability $cf $bf $eta $sb" >> "$JOBFILE"
                done
            done
        done
    done
done

TOTAL=$(wc -l < "$JOBFILE")
echo "Total jobs: $TOTAL (parallel=$PARALLEL, time_limit=${TLIMIT}s)"
echo "Estimated max wall time: $(( TOTAL * TLIMIT / PARALLEL / 3600 )) hours"
echo ""

# CSV header
echo "instance,branch,cut_freq,balas_freq,eta,sb,status,primal,dual,gap,nodes,lp_solves,time" > "$OUTFILE"

# Lock file for thread-safe CSV writes
LOCKFILE=$(mktemp)

run_one() {
    local f="$1" branch="$2" cf="$3" bf="$4" eta="$5" sb="$6"
    local name=$(basename "$f" .txt)
    local result
    result=$("$SOLVER" "$f" --time-limit "$TLIMIT" --verbosity 0 \
        --branch "$branch" --cut-frequency "$cf" --balas-frequency "$bf" \
        --reliability-eta "$eta" --reliability-sb "$sb" 2>&1)
    local status=$(echo "$result" | sed -n 's/^Status:[[:space:]]*//p' | tr -d '\r')
    local primal=$(echo "$result" | sed -n 's/^Primal:[[:space:]]*//p' | tr -d '\r')
    local dual=$(echo "$result" | sed -n 's/^Dual:[[:space:]]*//p' | tr -d '\r')
    local gap=$(echo "$result" | sed -n 's/^MIP gap:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')
    local nodes=$(echo "$result" | sed -n 's/^Nodes:[[:space:]]*//p' | tr -d '\r')
    local lp_solves=$(echo "$result" | sed -n 's/^LP solves:[[:space:]]*//p' | tr -d '\r')
    local wtime=$(echo "$result" | sed -n 's/^Wall time:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')
    local line="$name,$branch,$cf,$bf,$eta,$sb,$status,$primal,$dual,${gap:-0},$nodes,$lp_solves,$wtime"
    # Thread-safe write
    (
        flock 9
        echo "$line" >> "$OUTFILE"
    ) 9>"$LOCKFILE"
    echo "$line"
}

export -f run_one
export SOLVER TLIMIT OUTFILE LOCKFILE

DONE=0
while IFS= read -r job; do
    # Wait if we have $PARALLEL jobs running
    while [ "$(jobs -rp | wc -l)" -ge "$PARALLEL" ]; do
        wait -n 2>/dev/null || sleep 1
    done
    run_one $job &
    DONE=$((DONE + 1))
    # Progress every 10 jobs
    if [ $((DONE % 10)) -eq 0 ]; then
        echo "[Progress: $DONE / $TOTAL launched]"
    fi
done < "$JOBFILE"

# Wait for remaining jobs
wait

rm -f "$JOBFILE" "$LOCKFILE"

echo ""
echo "Done. $TOTAL jobs completed."
echo "Results written to $OUTFILE"
