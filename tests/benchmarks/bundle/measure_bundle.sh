#!/usr/bin/env bash
# =============================================================================
# Bundle Size Measurement for IQCP JCIM Manuscript
# =============================================================================
#
# Produces precise, reproducible bundle size measurements for the Performance
# Benchmarks section of the JCIM Application Note (Paper 2 / M2).
#
# Measures (raw and gzipped):
#   - Total dist/ directory size
#   - WASM module
#   - Initial JS entry point (main chunk)
#   - Lazy-loaded chunks (Three.js viewer3d, html2pdf, module-specific)
#   - Web Worker
#   - CSS
#   - Total transfer size (sum of all gzipped assets)
#
# Usage:
#   bash tests/benchmarks/bundle/measure_bundle.sh
#   bash tests/benchmarks/bundle/measure_bundle.sh --skip-build   # Use existing dist/
#
# Output:
#   - Human-readable table to stdout
#   - JSON to tests/benchmarks/bundle/results/bundle_sizes.json
# =============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
WEB_DIR="$REPO_ROOT/apps/web"
DIST_DIR="$WEB_DIR/dist"
ASSETS_DIR="$DIST_DIR/assets"
RESULTS_DIR="$REPO_ROOT/tests/benchmarks/bundle/results"
OUTPUT_JSON="$RESULTS_DIR/bundle_sizes.json"

SKIP_BUILD=false
if [[ "${1:-}" == "--skip-build" ]]; then
    SKIP_BUILD=true
fi

# Colors
if [ -t 1 ]; then
    BOLD='\033[1m'
    DIM='\033[2m'
    CYAN='\033[0;36m'
    GREEN='\033[0;32m'
    YELLOW='\033[0;33m'
    NC='\033[0m'
else
    BOLD=''
    DIM=''
    CYAN=''
    GREEN=''
    YELLOW=''
    NC=''
fi

# --- Helper functions --------------------------------------------------------

# Gzipped size in bytes for a single file
gz_bytes() {
    gzip -c "$1" | wc -c | tr -d ' '
}

# Raw size in bytes
raw_bytes() {
    wc -c < "$1" | tr -d ' '
}

# Format bytes as human-readable (KB with 1 decimal)
fmt_kb() {
    local bytes="$1"
    awk "BEGIN { printf \"%.1f\", $bytes / 1024 }"
}

# Format bytes as human-readable (MB with 2 decimals)
fmt_mb() {
    local bytes="$1"
    awk "BEGIN { printf \"%.2f\", $bytes / 1048576 }"
}

# Print a table row: name, raw_bytes, gz_bytes
table_row() {
    local name="$1"
    local raw="$2"
    local gz="$3"
    local raw_display gz_display
    if [ "$raw" -ge 1048576 ]; then
        raw_display="$(fmt_mb "$raw") MB"
    else
        raw_display="$(fmt_kb "$raw") KB"
    fi
    if [ "$gz" -ge 1048576 ]; then
        gz_display="$(fmt_mb "$gz") MB"
    else
        gz_display="$(fmt_kb "$gz") KB"
    fi
    printf "  %-35s %10s  %10s\n" "$name" "$raw_display" "$gz_display"
}

# --- Step 1: Build -----------------------------------------------------------

if [ "$SKIP_BUILD" = false ]; then
    echo ""
    printf "${BOLD}[1/4] Building production bundle...${NC}\n"
    echo ""

    # Check if WASM module exists
    if [ ! -f "$WEB_DIR/src/wasm/qc_wasm_bg.wasm" ]; then
        printf "${YELLOW}  WASM module not found. Building with wasm-pack...${NC}\n"
        (cd "$REPO_ROOT" && wasm-pack build crates/qc-wasm --release --target web --out-dir ../../apps/web/src/wasm 2>&1)
        echo ""
    fi

    (cd "$WEB_DIR" && npm run build 2>&1)
    echo ""
else
    echo ""
    printf "${BOLD}[1/4] Skipping build (--skip-build)${NC}\n"
fi

# Verify dist exists
if [ ! -d "$ASSETS_DIR" ]; then
    printf "\n${YELLOW}ERROR: $ASSETS_DIR not found. Run without --skip-build.${NC}\n"
    exit 1
fi

# --- Step 2: Collect metadata -------------------------------------------------

printf "${BOLD}[2/4] Collecting build metadata...${NC}\n"

DATE_ISO="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
DATE_SHORT="$(date +%Y-%m-%d)"
NODE_VERSION="$(node --version 2>/dev/null || echo 'unknown')"
VITE_VERSION="$(cd "$WEB_DIR" && npx vite --version 2>/dev/null | head -1 || echo 'unknown')"
WASM_PACK_VERSION="$(wasm-pack --version 2>/dev/null | awk '{print $2}' || echo 'unknown')"
RUST_VERSION="$(rustc --version 2>/dev/null | awk '{print $2}' || echo 'unknown')"
RUSTC_FULL="$(rustc --version 2>/dev/null || echo 'unknown')"

# --- Step 3: Measure sizes ----------------------------------------------------

printf "${BOLD}[3/4] Measuring bundle sizes...${NC}\n"
echo ""

# --- Identify files ---

# WASM module
WASM_FILE=""
for f in "$ASSETS_DIR"/qc_wasm_bg*.wasm; do
    [ -f "$f" ] && WASM_FILE="$f" && break
done

# Main JS chunk (largest index-*.js)
MAIN_JS=""
MAIN_JS_RAW=0
for f in "$ASSETS_DIR"/index-*.js; do
    [ -f "$f" ] || continue
    size=$(raw_bytes "$f")
    if [ "$size" -gt "$MAIN_JS_RAW" ]; then
        MAIN_JS="$f"
        MAIN_JS_RAW=$size
    fi
done

# Secondary index chunk (smaller index-*.js, e.g. dynamic import helper)
SECONDARY_JS=""
for f in "$ASSETS_DIR"/index-*.js; do
    [ -f "$f" ] || continue
    [ "$f" = "$MAIN_JS" ] && continue
    SECONDARY_JS="$f"
    break
done

# Three.js / viewer3d chunk
VIEWER3D_FILE=""
for f in "$ASSETS_DIR"/viewer3d-*.js; do
    [ -f "$f" ] && VIEWER3D_FILE="$f" && break
done

# html2pdf chunk
HTML2PDF_FILE=""
for f in "$ASSETS_DIR"/html2pdf-*.js; do
    [ -f "$f" ] && HTML2PDF_FILE="$f" && break
done

# Web Worker
WORKER_FILE=""
for f in "$ASSETS_DIR"/compute.worker-*.js; do
    [ -f "$f" ] && WORKER_FILE="$f" && break
done

# Module-specific lazy chunks (ModuleD, ModuleE, etc.)
MODULE_CHUNKS=()
for f in "$ASSETS_DIR"/Module*-*.js; do
    [ -f "$f" ] && MODULE_CHUNKS+=("$f")
done

# CSS
CSS_FILE=""
for f in "$ASSETS_DIR"/index-*.css; do
    [ -f "$f" ] && CSS_FILE="$f" && break
done

# --- Measure each file ---

# Accumulators for total transfer size
TOTAL_RAW=0
TOTAL_GZ=0

declare -A SIZES_RAW
declare -A SIZES_GZ

measure_file() {
    local key="$1"
    local file="$2"
    if [ -n "$file" ] && [ -f "$file" ]; then
        SIZES_RAW[$key]=$(raw_bytes "$file")
        SIZES_GZ[$key]=$(gz_bytes "$file")
        TOTAL_RAW=$(( TOTAL_RAW + ${SIZES_RAW[$key]} ))
        TOTAL_GZ=$(( TOTAL_GZ + ${SIZES_GZ[$key]} ))
    else
        SIZES_RAW[$key]=0
        SIZES_GZ[$key]=0
    fi
}

measure_file "wasm" "$WASM_FILE"
measure_file "initial_js" "$MAIN_JS"
measure_file "secondary_js" "$SECONDARY_JS"
measure_file "three_js_chunk" "$VIEWER3D_FILE"
measure_file "html2pdf_chunk" "$HTML2PDF_FILE"
measure_file "worker" "$WORKER_FILE"
measure_file "css" "$CSS_FILE"

# Module chunks
MODULE_CHUNK_ENTRIES=""
for mc in "${MODULE_CHUNKS[@]+"${MODULE_CHUNKS[@]}"}"; do
    [ -z "$mc" ] && continue
    mc_name="$(basename "$mc")"
    mc_raw=$(raw_bytes "$mc")
    mc_gz=$(gz_bytes "$mc")
    TOTAL_RAW=$(( TOTAL_RAW + mc_raw ))
    TOTAL_GZ=$(( TOTAL_GZ + mc_gz ))
    if [ -n "$MODULE_CHUNK_ENTRIES" ]; then
        MODULE_CHUNK_ENTRIES="$MODULE_CHUNK_ENTRIES,"
    fi
    MODULE_CHUNK_ENTRIES="$MODULE_CHUNK_ENTRIES
      {
        \"name\": \"$mc_name\",
        \"raw_bytes\": $mc_raw,
        \"gzipped_bytes\": $mc_gz
      }"
done

# Also measure font files (KaTeX, etc.) -- aggregate only
FONT_RAW=0
FONT_GZ=0
FONT_COUNT=0
for f in "$ASSETS_DIR"/*.woff2 "$ASSETS_DIR"/*.woff "$ASSETS_DIR"/*.ttf; do
    [ -f "$f" ] || continue
    fr=$(raw_bytes "$f")
    fg=$(gz_bytes "$f")
    FONT_RAW=$(( FONT_RAW + fr ))
    FONT_GZ=$(( FONT_GZ + fg ))
    FONT_COUNT=$(( FONT_COUNT + 1 ))
    TOTAL_RAW=$(( TOTAL_RAW + fr ))
    TOTAL_GZ=$(( TOTAL_GZ + fg ))
done

# Total dist directory (including index.html, headers, etc.)
DIST_TOTAL_RAW=0
while IFS= read -r -d '' f; do
    DIST_TOTAL_RAW=$(( DIST_TOTAL_RAW + $(raw_bytes "$f") ))
done < <(find "$DIST_DIR" -type f -print0)

# --- Print human-readable table ---

printf "${BOLD}=== IQCP Bundle Size Measurement ===${NC}\n"
printf "${DIM}  Date: %s | Node: %s | Rust: %s${NC}\n" "$DATE_SHORT" "$NODE_VERSION" "$RUST_VERSION"
echo ""

printf "  %-35s %10s  %10s\n" "ASSET" "RAW" "GZIPPED"
printf "  %-35s %10s  %10s\n" "-----------------------------------" "----------" "----------"

# Core assets
table_row "Initial JS (main chunk)" "${SIZES_RAW[initial_js]}" "${SIZES_GZ[initial_js]}"
if [ "${SIZES_RAW[secondary_js]}" -gt 0 ]; then
    table_row "Secondary JS (index helper)" "${SIZES_RAW[secondary_js]}" "${SIZES_GZ[secondary_js]}"
fi
table_row "WASM module (qc_wasm_bg)" "${SIZES_RAW[wasm]}" "${SIZES_GZ[wasm]}"
table_row "CSS" "${SIZES_RAW[css]}" "${SIZES_GZ[css]}"
table_row "Web Worker" "${SIZES_RAW[worker]}" "${SIZES_GZ[worker]}"

echo ""
printf "  ${DIM}Lazy-loaded chunks:${NC}\n"
table_row "  Three.js chunk (viewer3d)" "${SIZES_RAW[three_js_chunk]}" "${SIZES_GZ[three_js_chunk]}"
table_row "  html2pdf" "${SIZES_RAW[html2pdf_chunk]}" "${SIZES_GZ[html2pdf_chunk]}"

for mc in "${MODULE_CHUNKS[@]+"${MODULE_CHUNKS[@]}"}"; do
    [ -z "$mc" ] && continue
    mc_name="$(basename "$mc")"
    mc_raw=$(raw_bytes "$mc")
    mc_gz=$(gz_bytes "$mc")
    table_row "  $mc_name" "$mc_raw" "$mc_gz"
done

if [ "$FONT_COUNT" -gt 0 ]; then
    echo ""
    printf "  ${DIM}Fonts (%d files):${NC}\n" "$FONT_COUNT"
    table_row "  KaTeX fonts (aggregate)" "$FONT_RAW" "$FONT_GZ"
fi

echo ""
printf "  %-35s %10s  %10s\n" "-----------------------------------" "----------" "----------"
table_row "TOTAL (all assets)" "$TOTAL_RAW" "$TOTAL_GZ"

# dist/ total is raw-only (no aggregate gzip)
local_dist_display=$(fmt_mb "$DIST_TOTAL_RAW")
printf "  %-35s %7s MB  %10s\n" "TOTAL dist/ directory" "$local_dist_display" "n/a"

echo ""

# --- Budget summary ---
printf "${BOLD}Budget Comparison (gzipped):${NC}\n"

check_budget() {
    local name="$1"
    local actual_gz="$2"
    local budget_kb="$3"
    local budget_bytes=$(( budget_kb * 1024 ))
    local actual_kb
    actual_kb=$(awk "BEGIN { printf \"%.1f\", $actual_gz / 1024 }")

    if [ "$actual_gz" -le "$budget_bytes" ]; then
        local headroom
        headroom=$(awk "BEGIN { printf \"%.1f\", ($budget_bytes - $actual_gz) / 1024 }")
        printf "  ${GREEN}PASS${NC}  %-30s %8s KB  (budget: %d KB, headroom: %s KB)\n" "$name" "$actual_kb" "$budget_kb" "$headroom"
    else
        local over
        over=$(awk "BEGIN { printf \"%.1f\", ($actual_gz - $budget_bytes) / 1024 }")
        printf "  ${YELLOW}OVER${NC}  %-30s %8s KB  (budget: %d KB, over: %s KB)\n" "$name" "$actual_kb" "$budget_kb" "$over"
    fi
}

check_budget "Initial JS" "${SIZES_GZ[initial_js]}" 600
check_budget "WASM module" "${SIZES_GZ[wasm]}" 250
check_budget "Three.js chunk" "${SIZES_GZ[three_js_chunk]}" 300

echo ""

# --- Step 4: Write JSON output ------------------------------------------------

printf "${BOLD}[4/4] Writing results to $OUTPUT_JSON${NC}\n"

# Build chunks array
CHUNKS_JSON="["
CHUNK_FIRST=true

add_chunk() {
    local name="$1"
    local raw="$2"
    local gz="$3"
    local category="$4"
    if [ "$CHUNK_FIRST" = true ]; then
        CHUNK_FIRST=false
    else
        CHUNKS_JSON="$CHUNKS_JSON,"
    fi
    CHUNKS_JSON="$CHUNKS_JSON
    {
      \"name\": \"$name\",
      \"category\": \"$category\",
      \"raw_bytes\": $raw,
      \"gzipped_bytes\": $gz
    }"
}

# Add all chunks
if [ -n "$MAIN_JS" ]; then
    add_chunk "$(basename "$MAIN_JS")" "${SIZES_RAW[initial_js]}" "${SIZES_GZ[initial_js]}" "initial"
fi
if [ "${SIZES_RAW[secondary_js]}" -gt 0 ]; then
    add_chunk "$(basename "$SECONDARY_JS")" "${SIZES_RAW[secondary_js]}" "${SIZES_GZ[secondary_js]}" "initial"
fi
if [ -n "$WASM_FILE" ]; then
    add_chunk "$(basename "$WASM_FILE")" "${SIZES_RAW[wasm]}" "${SIZES_GZ[wasm]}" "wasm"
fi
if [ -n "$CSS_FILE" ]; then
    add_chunk "$(basename "$CSS_FILE")" "${SIZES_RAW[css]}" "${SIZES_GZ[css]}" "css"
fi
if [ -n "$WORKER_FILE" ]; then
    add_chunk "$(basename "$WORKER_FILE")" "${SIZES_RAW[worker]}" "${SIZES_GZ[worker]}" "worker"
fi
if [ -n "$VIEWER3D_FILE" ]; then
    add_chunk "$(basename "$VIEWER3D_FILE")" "${SIZES_RAW[three_js_chunk]}" "${SIZES_GZ[three_js_chunk]}" "lazy"
fi
if [ -n "$HTML2PDF_FILE" ]; then
    add_chunk "$(basename "$HTML2PDF_FILE")" "${SIZES_RAW[html2pdf_chunk]}" "${SIZES_GZ[html2pdf_chunk]}" "lazy"
fi
for mc in "${MODULE_CHUNKS[@]+"${MODULE_CHUNKS[@]}"}"; do
    [ -z "$mc" ] && continue
    mc_name="$(basename "$mc")"
    mc_raw=$(raw_bytes "$mc")
    mc_gz=$(gz_bytes "$mc")
    add_chunk "$mc_name" "$mc_raw" "$mc_gz" "lazy"
done

CHUNKS_JSON="$CHUNKS_JSON
  ]"

# Write final JSON
cat > "$OUTPUT_JSON" << ENDJSON
{
  "metadata": {
    "date": "$DATE_ISO",
    "date_short": "$DATE_SHORT",
    "node_version": "$NODE_VERSION",
    "vite_version": "$VITE_VERSION",
    "wasm_pack_version": "$WASM_PACK_VERSION",
    "rust_version": "$RUST_VERSION",
    "rustc": "$RUSTC_FULL",
    "measurement_tool": "tests/benchmarks/bundle/measure_bundle.sh",
    "compression": "gzip (default level)"
  },
  "sizes": {
    "wasm": {
      "raw_bytes": ${SIZES_RAW[wasm]},
      "gzipped_bytes": ${SIZES_GZ[wasm]}
    },
    "initial_js": {
      "raw_bytes": ${SIZES_RAW[initial_js]},
      "gzipped_bytes": ${SIZES_GZ[initial_js]}
    },
    "three_js_chunk": {
      "raw_bytes": ${SIZES_RAW[three_js_chunk]},
      "gzipped_bytes": ${SIZES_GZ[three_js_chunk]}
    },
    "css": {
      "raw_bytes": ${SIZES_RAW[css]},
      "gzipped_bytes": ${SIZES_GZ[css]}
    },
    "worker": {
      "raw_bytes": ${SIZES_RAW[worker]},
      "gzipped_bytes": ${SIZES_GZ[worker]}
    },
    "html2pdf_chunk": {
      "raw_bytes": ${SIZES_RAW[html2pdf_chunk]},
      "gzipped_bytes": ${SIZES_GZ[html2pdf_chunk]}
    },
    "fonts": {
      "count": $FONT_COUNT,
      "raw_bytes": $FONT_RAW,
      "gzipped_bytes": $FONT_GZ
    },
    "total_transfer": {
      "raw_bytes": $TOTAL_RAW,
      "gzipped_bytes": $TOTAL_GZ
    },
    "total_dist": {
      "raw_bytes": $DIST_TOTAL_RAW
    }
  },
  "budgets": {
    "initial_js_kb": 600,
    "wasm_kb": 250,
    "three_js_chunk_kb": 300
  },
  "chunks": $CHUNKS_JSON
}
ENDJSON

printf "  Saved: ${CYAN}%s${NC}\n" "$OUTPUT_JSON"
echo ""
printf "${GREEN}${BOLD}Measurement complete.${NC}\n"
echo ""
