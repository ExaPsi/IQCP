#!/usr/bin/env bash
#
# Run all PySCF benchmark scripts for JCIM M2 manuscript.
#
# Usage:
#   bash tests/benchmarks/pyscf/run_all.sh
#
# Requires the PySCF virtual environment at .venv/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
VENV="$PROJECT_ROOT/.venv"

echo "========================================"
echo "  PySCF Benchmark Suite"
echo "========================================"
echo ""

# Activate virtual environment
if [ ! -d "$VENV" ]; then
    echo "ERROR: Virtual environment not found at $VENV"
    echo "Create it:  python -m venv .venv && pip install pyscf"
    exit 1
fi

source "$VENV/bin/activate"
echo "Python: $(python --version)"
echo "PySCF:  $(python -c 'import pyscf; print(pyscf.__version__)')"
echo ""

# Run SCF benchmarks
echo "========================================="
echo "  1/3  SCF Timing Benchmarks"
echo "========================================="
python "$SCRIPT_DIR/benchmark_scf.py"
echo ""

# Run gradient benchmarks
echo "========================================="
echo "  2/3  Gradient Timing Benchmarks"
echo "========================================="
python "$SCRIPT_DIR/benchmark_gradients.py"
echo ""

# Run optimization benchmarks
echo "========================================="
echo "  3/3  Optimization Timing Benchmarks"
echo "========================================="
python "$SCRIPT_DIR/benchmark_optimization.py"
echo ""

echo "========================================="
echo "  All benchmarks complete."
echo "  Results in: $SCRIPT_DIR/results/"
echo "========================================="
ls -la "$SCRIPT_DIR/results/"
