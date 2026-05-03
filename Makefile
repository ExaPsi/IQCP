# =============================================================================
# IQCP — Reviewer-friendly build targets
# =============================================================================
#
# Quick start for reviewers:
#   make pyscf-env   # one-time: create .venv with PySCF + dependencies
#   make all         # check + test + wasm + web
#
# Individual targets:
#   make check       # cargo fmt --check + clippy
#   make test        # cargo test --workspace
#   make bench-build # compile benches without running
#   make wasm        # build qc-wasm and qc-wasm-spectra modules
#   make web         # install npm deps and build the SPA
#   make clean       # remove all build artifacts
# =============================================================================

.PHONY: help all check test bench-build wasm wasm-core wasm-spectra web pyscf-env clean

WEB_DIR := apps/web
WASM_OUT := $(WEB_DIR)/src/wasm
VENV := .venv

# Default target: print help
help:
	@echo "IQCP build targets:"
	@echo "  make check        — cargo fmt --check + cargo clippy --workspace -- -D warnings"
	@echo "  make test         — cargo test --workspace"
	@echo "  make bench-build  — cargo bench --no-run (compile benches without running)"
	@echo "  make wasm         — wasm-pack build for qc-wasm and qc-wasm-spectra"
	@echo "  make web          — npm install + npm run build in apps/web"
	@echo "  make pyscf-env    — create .venv and install pinned PySCF dependencies"
	@echo "  make all          — check + test + wasm + web"
	@echo "  make clean        — cargo clean + remove web build artifacts"

# -----------------------------------------------------------------------------
# Rust checks
# -----------------------------------------------------------------------------
check:
	@echo ">>> Checking Rust formatting (cargo fmt --all -- --check)"
	cargo fmt --all -- --check
	@echo ">>> Linting Rust workspace (cargo clippy --workspace -- -D warnings)"
	cargo clippy --workspace -- -D warnings

test:
	@echo ">>> Running Rust workspace tests (cargo test --workspace)"
	cargo test --workspace

bench-build:
	@echo ">>> Compiling benchmarks without running (cargo bench --no-run)"
	cargo bench --workspace --no-run

# -----------------------------------------------------------------------------
# WASM build (requires wasm-pack: https://rustwasm.github.io/wasm-pack/)
# -----------------------------------------------------------------------------
wasm: wasm-core wasm-spectra

wasm-core:
	@echo ">>> Building qc-wasm module (release, target=web)"
	wasm-pack build crates/qc-wasm --release --target web --out-dir ../../$(WASM_OUT)

wasm-spectra:
	@echo ">>> Building qc-wasm-spectra module (release, target=web)"
	wasm-pack build crates/qc-wasm-spectra --release --target web --out-dir ../../$(WASM_OUT)/spectra

# -----------------------------------------------------------------------------
# Web frontend
# -----------------------------------------------------------------------------
web:
	@echo ">>> Installing web dependencies (npm install in $(WEB_DIR))"
	cd $(WEB_DIR) && npm install
	@echo ">>> Building production bundle (npm run build in $(WEB_DIR))"
	cd $(WEB_DIR) && npm run build

# -----------------------------------------------------------------------------
# PySCF validation environment
# -----------------------------------------------------------------------------
pyscf-env:
	@echo ">>> Creating Python virtual environment at $(VENV)"
	python3 -m venv $(VENV)
	@echo ">>> Installing pinned PySCF dependencies from requirements.txt"
	$(VENV)/bin/pip install --upgrade pip
	$(VENV)/bin/pip install -r requirements.txt
	@echo ">>> PySCF environment ready. Activate with: source $(VENV)/bin/activate"

# -----------------------------------------------------------------------------
# Aggregate target
# -----------------------------------------------------------------------------
all: check test wasm web
	@echo ">>> All build targets completed successfully."

# -----------------------------------------------------------------------------
# Cleanup
# -----------------------------------------------------------------------------
clean:
	@echo ">>> Cleaning Rust build artifacts (cargo clean)"
	cargo clean
	@echo ">>> Removing web build outputs ($(WEB_DIR)/dist, $(WEB_DIR)/node_modules/.vite)"
	rm -rf $(WEB_DIR)/dist $(WEB_DIR)/node_modules/.vite
