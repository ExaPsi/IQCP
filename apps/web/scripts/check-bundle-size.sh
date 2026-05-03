#!/usr/bin/env bash
# =============================================================================
# Bundle Size Budget Check for IQCP Phase 2
# =============================================================================
#
# Enforces bundle size budgets and code-split verification.
# Run after: npm run build (or use npm run check:bundle which builds first)
#
# Usage:
#   bash scripts/check-bundle-size.sh          # Check existing build
#   bash scripts/check-bundle-size.sh --build  # Build first, then check
#
# Budget limits (gzipped):
#   Initial JS:     600 KB (WARNING if exceeded -- Plotly.js known issue)
#   WASM module:    250 KB (HARD FAIL if exceeded)
#   Three.js chunk: 300 KB (HARD FAIL if exceeded)
#
# Code-split rule:
#   Main JS chunk must NOT contain Three.js / R3F / Drei code.
#
# Exit codes:
#   0  All hard budgets pass and code-split verified
#   1  At least one hard budget exceeded or code-split violated
# =============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------

DIST_DIR="dist/assets"
JS_BUDGET_KB=600
WASM_BUDGET_KB=250
THREEJS_BUDGET_KB=300

# Counters
PASS=0
FAIL=0
WARN=0

# Colors (if terminal supports them)
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    CYAN='\033[0;36m'
    BOLD='\033[1m'
    NC='\033[0m'
else
    GREEN=''
    RED=''
    YELLOW=''
    CYAN=''
    BOLD=''
    NC=''
fi

# --- Helper functions --------------------------------------------------------

gz_size_bytes() {
    gzip -c "$1" | wc -c
}

gz_size_kb() {
    local bytes
    bytes=$(gz_size_bytes "$1")
    echo $(( bytes / 1024 ))
}

report_pass() {
    printf "${GREEN}  PASS${NC}  %-30s %6s KB gz  (budget: %s KB)\n" "$1" "$2" "$3"
    PASS=$(( PASS + 1 ))
}

report_fail() {
    printf "${RED}  FAIL${NC}  %-30s %6s KB gz  (budget: %s KB) ${RED}OVER by %s KB${NC}\n" "$1" "$2" "$3" "$4"
    FAIL=$(( FAIL + 1 ))
}

report_warn() {
    printf "${YELLOW}  WARN${NC}  %-30s %6s KB gz  (budget: %s KB) ${YELLOW}OVER by %s KB${NC}\n" "$1" "$2" "$3" "$4"
    WARN=$(( WARN + 1 ))
}

report_info() {
    printf "${CYAN}  INFO${NC}  %-30s %6s KB gz\n" "$1" "$2"
}

# --- Optionally build first --------------------------------------------------

if [[ "${1:-}" == "--build" ]]; then
    echo ""
    printf "${BOLD}Building production bundle...${NC}\n"
    npm run build
    echo ""
fi

# --- Verify dist directory exists ---------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DIST_PATH="$PROJECT_DIR/$DIST_DIR"

if [ ! -d "$DIST_PATH" ]; then
    printf "${RED}ERROR: dist/assets/ not found. Run 'npm run build' first.${NC}\n"
    exit 1
fi

# --- Find chunks --------------------------------------------------------------

echo ""
printf "${BOLD}=== IQCP Bundle Size Budget Check ===${NC}\n"
echo ""

# Find main JS chunk (largest index-*.js file)
MAIN_JS=""
MAIN_JS_SIZE=0
for f in "$DIST_PATH"/index-*.js; do
    [ -f "$f" ] || continue
    size=$(wc -c < "$f")
    if [ "$size" -gt "$MAIN_JS_SIZE" ]; then
        MAIN_JS="$f"
        MAIN_JS_SIZE=$size
    fi
done

# Find WASM module
WASM_FILE=""
for f in "$DIST_PATH"/qc_wasm_bg-*.wasm "$DIST_PATH"/qc_wasm_bg.wasm; do
    [ -f "$f" ] && WASM_FILE="$f" && break
done

# Find Three.js / viewer3d chunk
VIEWER3D_FILE=""
for f in "$DIST_PATH"/viewer3d-*.js; do
    [ -f "$f" ] && VIEWER3D_FILE="$f" && break
done

# Find html2pdf chunk (lazy-loaded, informational only)
HTML2PDF_FILE=""
for f in "$DIST_PATH"/html2pdf-*.js; do
    [ -f "$f" ] && HTML2PDF_FILE="$f" && break
done

# Find worker chunk (informational only)
WORKER_FILE=""
for f in "$DIST_PATH"/compute.worker-*.js; do
    [ -f "$f" ] && WORKER_FILE="$f" && break
done

# --- Check budgets ------------------------------------------------------------

printf "${BOLD}Budget Checks:${NC}\n"

# 1. Main JS bundle
if [ -n "$MAIN_JS" ]; then
    MAIN_GZ_KB=$(gz_size_kb "$MAIN_JS")
    if [ "$MAIN_GZ_KB" -le "$JS_BUDGET_KB" ]; then
        report_pass "Main JS (initial)" "$MAIN_GZ_KB" "$JS_BUDGET_KB"
    else
        OVER=$(( MAIN_GZ_KB - JS_BUDGET_KB ))
        # WARNING only -- Plotly.js is a known contributor that needs a separate
        # migration story (potential uPlot replacement).
        report_warn "Main JS (initial)" "$MAIN_GZ_KB" "$JS_BUDGET_KB" "$OVER"
        printf "         ${YELLOW}Note: Plotly.js dominates main chunk. Track uPlot migration separately.${NC}\n"
    fi
else
    printf "${RED}  FAIL${NC}  Main JS chunk not found in $DIST_DIR\n"
    FAIL=$(( FAIL + 1 ))
fi

# 2. WASM module
if [ -n "$WASM_FILE" ]; then
    WASM_GZ_KB=$(gz_size_kb "$WASM_FILE")
    if [ "$WASM_GZ_KB" -le "$WASM_BUDGET_KB" ]; then
        report_pass "WASM module" "$WASM_GZ_KB" "$WASM_BUDGET_KB"
    else
        OVER=$(( WASM_GZ_KB - WASM_BUDGET_KB ))
        report_fail "WASM module" "$WASM_GZ_KB" "$WASM_BUDGET_KB" "$OVER"
    fi
else
    printf "${RED}  FAIL${NC}  WASM module not found in $DIST_DIR\n"
    FAIL=$(( FAIL + 1 ))
fi

# 3. Three.js / viewer3d chunk
if [ -n "$VIEWER3D_FILE" ]; then
    V3D_GZ_KB=$(gz_size_kb "$VIEWER3D_FILE")
    if [ "$V3D_GZ_KB" -le "$THREEJS_BUDGET_KB" ]; then
        report_pass "Three.js chunk (viewer3d)" "$V3D_GZ_KB" "$THREEJS_BUDGET_KB"
    else
        OVER=$(( V3D_GZ_KB - THREEJS_BUDGET_KB ))
        report_fail "Three.js chunk (viewer3d)" "$V3D_GZ_KB" "$THREEJS_BUDGET_KB" "$OVER"
    fi
else
    # viewer3d chunk is optional -- only exists if Three.js deps are installed
    printf "${CYAN}  INFO${NC}  viewer3d chunk not found (Three.js not bundled)\n"
fi

# --- Informational: other chunks ----------------------------------------------

echo ""
printf "${BOLD}Other Chunks (informational):${NC}\n"

if [ -n "$HTML2PDF_FILE" ]; then
    report_info "html2pdf (lazy)" "$(gz_size_kb "$HTML2PDF_FILE")"
fi

if [ -n "$WORKER_FILE" ]; then
    report_info "Web Worker" "$(gz_size_kb "$WORKER_FILE")"
fi

# WASM JS glue
for f in "$DIST_PATH"/qc_wasm-*.js; do
    [ -f "$f" ] && report_info "WASM glue JS" "$(gz_size_kb "$f")" && break
done

# Worker helpers
for f in "$DIST_PATH"/workerHelpers-*.js; do
    [ -f "$f" ] && report_info "Worker helpers" "$(gz_size_kb "$f")" && break
done

# --- Code-split verification --------------------------------------------------

echo ""
printf "${BOLD}Code-Split Verification:${NC}\n"

if [ -n "$MAIN_JS" ]; then
    SPLIT_OK=true

    # Check for Three.js markers in main chunk
    if grep -q 'WebGLRenderer' "$MAIN_JS" 2>/dev/null; then
        printf "${RED}  FAIL${NC}  Main chunk contains 'WebGLRenderer' (Three.js leak)\n"
        SPLIT_OK=false
        FAIL=$(( FAIL + 1 ))
    fi

    if grep -q '@react-three' "$MAIN_JS" 2>/dev/null; then
        printf "${RED}  FAIL${NC}  Main chunk contains '@react-three' (R3F/Drei leak)\n"
        SPLIT_OK=false
        FAIL=$(( FAIL + 1 ))
    fi

    # Check for Three.js module-level exports that would indicate bundling
    if grep -q 'THREE\.Scene\|THREE\.Mesh\|THREE\.Vector3' "$MAIN_JS" 2>/dev/null; then
        printf "${RED}  FAIL${NC}  Main chunk contains THREE.* namespace references\n"
        SPLIT_OK=false
        FAIL=$(( FAIL + 1 ))
    fi

    if $SPLIT_OK; then
        printf "${GREEN}  PASS${NC}  Three.js / R3F / Drei not found in main chunk\n"
        PASS=$(( PASS + 1 ))
    fi
else
    printf "${YELLOW}  SKIP${NC}  Cannot verify code-split (main chunk not found)\n"
fi

# --- Summary ------------------------------------------------------------------

echo ""
printf "${BOLD}=== Summary ===${NC}\n"
TOTAL=$(( PASS + FAIL + WARN ))
printf "  Total checks: %d\n" "$TOTAL"
printf "  ${GREEN}PASS: %d${NC}  ${RED}FAIL: %d${NC}  ${YELLOW}WARN: %d${NC}\n" "$PASS" "$FAIL" "$WARN"

if [ "$FAIL" -gt 0 ]; then
    echo ""
    printf "${RED}${BOLD}Bundle budget check FAILED.${NC}\n"
    printf "Fix the failing checks before committing.\n"
    exit 1
else
    echo ""
    if [ "$WARN" -gt 0 ]; then
        printf "${YELLOW}${BOLD}Bundle budget check PASSED with warnings.${NC}\n"
        printf "Warnings do not block commits but should be addressed.\n"
    else
        printf "${GREEN}${BOLD}Bundle budget check PASSED.${NC}\n"
    fi
    exit 0
fi
