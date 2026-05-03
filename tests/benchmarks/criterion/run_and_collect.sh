#!/usr/bin/env bash
# run_and_collect.sh — Run Criterion benchmarks and collect results into JSON
#
# Usage:
#   bash tests/benchmarks/criterion/run_and_collect.sh
#
# Output:
#   tests/benchmarks/criterion/results_native.json
#
# This script:
#   1. Runs all qc-core Criterion benchmarks (cargo bench)
#   2. Extracts timing data from Criterion's JSON output
#   3. Collects system information (CPU, Rust version, date)
#   4. Writes consolidated results to results_native.json

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
OUTPUT_FILE="$SCRIPT_DIR/results_native.json"
CRITERION_DIR="$PROJECT_ROOT/target/criterion"

echo "=== IQCP Criterion Benchmark Runner ==="
echo "Project root: $PROJECT_ROOT"
echo "Output file:  $OUTPUT_FILE"
echo ""

# 1. Collect system information
echo "--- Collecting system info ---"
RUST_VERSION=$(rustc --version 2>/dev/null || echo "unknown")
CARGO_VERSION=$(cargo --version 2>/dev/null || echo "unknown")
DATE_ISO=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
HOSTNAME=$(hostname 2>/dev/null || echo "unknown")
OS_INFO=$(uname -srm 2>/dev/null || echo "unknown")

# CPU info (Linux)
if [ -f /proc/cpuinfo ]; then
    CPU_MODEL=$(grep -m1 "model name" /proc/cpuinfo 2>/dev/null | cut -d: -f2 | xargs || echo "unknown")
    CPU_CORES=$(grep -c "^processor" /proc/cpuinfo 2>/dev/null || echo "unknown")
else
    CPU_MODEL=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "unknown")
    CPU_CORES=$(sysctl -n hw.ncpu 2>/dev/null || echo "unknown")
fi

echo "  Rust:  $RUST_VERSION"
echo "  CPU:   $CPU_MODEL ($CPU_CORES cores)"
echo "  OS:    $OS_INFO"
echo ""

# 2. Run benchmarks
echo "--- Running benchmarks ---"
cd "$PROJECT_ROOT"
cargo bench --package qc-core 2>&1 | tee /tmp/criterion_bench_output.txt
echo ""
echo "--- Benchmarks complete ---"
echo ""

# 3. Extract results from Criterion JSON output
echo "--- Extracting results ---"

# Start building the JSON output
{
    echo '{'
    echo '  "metadata": {'
    echo "    \"date\": \"$DATE_ISO\","
    echo "    \"rust_version\": \"$RUST_VERSION\","
    echo "    \"cargo_version\": \"$CARGO_VERSION\","
    echo "    \"hostname\": \"$HOSTNAME\","
    echo "    \"os\": \"$OS_INFO\","
    echo "    \"cpu_model\": \"$CPU_MODEL\","
    echo "    \"cpu_cores\": \"$CPU_CORES\""
    echo '  },'
    echo '  "benchmarks": ['

    first=true

    # Iterate over Criterion benchmark directories
    if [ -d "$CRITERION_DIR" ]; then
        for estimates_file in $(find "$CRITERION_DIR" -name "estimates.json" -path "*/new/*" 2>/dev/null | sort); do
            # Extract benchmark name from directory structure
            # Path: target/criterion/<group>/<bench_name>/new/estimates.json
            # or:   target/criterion/<bench_name>/new/estimates.json
            bench_path="${estimates_file#$CRITERION_DIR/}"
            bench_name="${bench_path%/new/estimates.json}"

            # Parse timing from estimates.json (values are in nanoseconds)
            if [ -f "$estimates_file" ]; then
                mean_ns=$(python3 -c "
import json, sys
with open('$estimates_file') as f:
    d = json.load(f)
print(d.get('mean', {}).get('point_estimate', 0))
" 2>/dev/null || echo "0")

                std_dev_ns=$(python3 -c "
import json, sys
with open('$estimates_file') as f:
    d = json.load(f)
print(d.get('std_dev', {}).get('point_estimate', 0))
" 2>/dev/null || echo "0")

                median_ns=$(python3 -c "
import json, sys
with open('$estimates_file') as f:
    d = json.load(f)
print(d.get('median', {}).get('point_estimate', 0))
" 2>/dev/null || echo "0")

                if [ "$first" = true ]; then
                    first=false
                else
                    echo '    ,'
                fi

                echo '    {'
                echo "      \"name\": \"$bench_name\","
                echo "      \"mean_ns\": $mean_ns,"
                echo "      \"std_dev_ns\": $std_dev_ns,"
                echo "      \"median_ns\": $median_ns"
                echo -n '    }'
            fi
        done
    fi

    echo ''
    echo '  ]'
    echo '}'
} > "$OUTPUT_FILE"

echo "Results written to: $OUTPUT_FILE"
echo ""

# 4. Print summary
echo "=== Benchmark Summary ==="
if command -v python3 &>/dev/null && [ -f "$OUTPUT_FILE" ]; then
    python3 -c "
import json

with open('$OUTPUT_FILE') as f:
    data = json.load(f)

print(f\"Date: {data['metadata']['date']}\")
print(f\"Rust: {data['metadata']['rust_version']}\")
print(f\"CPU:  {data['metadata']['cpu_model']}\")
print()
print(f\"{'Benchmark':<50} {'Mean':>12} {'Std Dev':>12}\")
print('-' * 76)

for b in data['benchmarks']:
    mean = b['mean_ns']
    std = b['std_dev_ns']

    # Format with appropriate units
    if mean >= 1e9:
        mean_str = f\"{mean/1e9:.3f} s\"
        std_str = f\"{std/1e9:.3f} s\"
    elif mean >= 1e6:
        mean_str = f\"{mean/1e6:.3f} ms\"
        std_str = f\"{std/1e6:.3f} ms\"
    elif mean >= 1e3:
        mean_str = f\"{mean/1e3:.3f} us\"
        std_str = f\"{std/1e3:.3f} us\"
    else:
        mean_str = f\"{mean:.1f} ns\"
        std_str = f\"{std:.1f} ns\"

    print(f\"{b['name']:<50} {mean_str:>12} {std_str:>12}\")
" 2>/dev/null || echo "(install python3 for formatted summary)"
fi

echo ""
echo "=== Done ==="
