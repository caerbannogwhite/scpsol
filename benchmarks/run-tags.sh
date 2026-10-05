#!/bin/bash
# Benchmark all tagged solver versions across SCP instance sets.
#
# Usage: bash benchmarks/run-tags.sh [time_limit]
#   time_limit: per-instance time limit in seconds (default: 300)
#
# Outputs:
#   benchmarks/results/<tag>.csv   — per-instance results for each tag
#   benchmarks/results/summary.csv — combined summary table

set -e

TLIMIT="${1:-300}"
RESULTS_DIR="benchmarks/results"
mkdir -p "$RESULTS_DIR"

TAGS=(v0.1-baseline v0.2-cuts-balas v0.3-reliability v0.4-preprocess v0.5-fast-preprocess v0.6-best-bound v0.7-dominance-finder v0.8-current)

# Collect instance files
INSTANCES=""
for pat in data/scp4{1,2,3,4,5,6,7,8,9,10}.txt data/scp5{1,2,3,4,5,6,7,8,9,10}.txt data/scpa{1,2,3,4,5}.txt data/scpb{1,2,3,4,5}.txt data/scpnre{1,2,3,4,5}.txt data/scpnrf{1,2,3,4,5}.txt; do
    if [ -f "$pat" ]; then
        INSTANCES="$INSTANCES $pat"
    fi
done

echo "=== Cross-version SCP benchmark ==="
echo "Tags: ${TAGS[*]}"
echo "Time limit: ${TLIMIT}s per instance"
echo "Instances: $(echo $INSTANCES | wc -w) files"
echo ""

# Save current branch/commit to restore later
ORIG_REF=$(git rev-parse HEAD)
ORIG_BRANCH=$(git branch --show-current 2>/dev/null || echo "")

for tag in "${TAGS[@]}"; do
    echo "========================================"
    echo "Building $tag ..."
    echo "========================================"

    # Checkout tag
    git checkout "$tag" --quiet 2>/dev/null

    # Build in a tag-specific build directory
    BUILD_DIR="build-${tag}"
    cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DFAST_BUILD=ON 2>/dev/null | tail -1
    cmake --build "$BUILD_DIR" --config Release 2>&1 | tail -3

    SOLVER="${BUILD_DIR}/Release/scpsol.exe"
    if [ ! -f "$SOLVER" ]; then
        echo "ERROR: $SOLVER not found, skipping $tag"
        continue
    fi

    OUTFILE="${RESULTS_DIR}/${tag}.csv"
    echo "instance,status,primal,dual,mip_gap_pct,nodes,lp_solves,wall_time_s" > "$OUTFILE"

    echo ""
    echo "Running $tag on all instances (limit=${TLIMIT}s)..."
    for f in $INSTANCES; do
        name=$(basename "$f" .txt)
        printf "  %-10s " "$name"

        # Use absolute path for data file since we may be at different checkout
        ABS_F="$(cd "$(dirname "$f")" && pwd)/$(basename "$f")"

        # Retry runs whose output is empty or incomplete (see issue #2)
        result=""
        for attempt in 1 2 3; do
            result=$("$SOLVER" "$ABS_F" --time-limit "$TLIMIT" --verbosity 0 2>&1) || true
            if echo "$result" | grep -q '^Status:'; then
                break
            fi
            result=""
        done
        status=$(echo "$result" | sed -n 's/^Status:[[:space:]]*//p' | tr -d '\r')
        primal=$(echo "$result" | sed -n 's/^Primal:[[:space:]]*//p' | tr -d '\r')
        dual=$(echo "$result" | sed -n 's/^Dual:[[:space:]]*//p' | tr -d '\r')
        gap=$(echo "$result" | sed -n 's/^MIP gap:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')
        nodes=$(echo "$result" | sed -n 's/^Nodes:[[:space:]]*//p' | tr -d '\r')
        lp_solves=$(echo "$result" | sed -n 's/^LP solves:[[:space:]]*//p' | tr -d '\r')
        wtime=$(echo "$result" | sed -n 's/^Wall time:[[:space:]]*\([0-9.]*\).*/\1/p' | tr -d '\r')

        # Handle missing fields for older solver versions
        status="${status:-Unknown}"
        primal="${primal:-inf}"
        dual="${dual:-inf}"
        gap="${gap:-inf}"
        nodes="${nodes:-0}"
        lp_solves="${lp_solves:-0}"
        wtime="${wtime:-${TLIMIT}}"

        printf "status=%-10s primal=%-6s gap=%-10s nodes=%-6s time=%s s\n" \
            "$status" "$primal" "${gap}%" "$nodes" "$wtime"
        echo "$name,$status,$primal,$dual,$gap,$nodes,$lp_solves,$wtime" >> "$OUTFILE"
    done

    echo ""
    echo "Results written to $OUTFILE"
    echo ""
done

# Restore original branch
if [ -n "$ORIG_BRANCH" ]; then
    git checkout "$ORIG_BRANCH" --quiet 2>/dev/null
else
    git checkout "$ORIG_REF" --quiet 2>/dev/null
fi

echo "========================================"
echo "All benchmarks complete."
echo "Results in $RESULTS_DIR/"
echo "========================================"
