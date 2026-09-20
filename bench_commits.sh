#!/usr/bin/env bash
set -e

# The optimization sequence to benchmark:
COMMITS=(
    "abr-v0.4.0"     # Historical v0.4.0 release
    "main"           # Branch base (c912ce51)
    "v1.0-release"   # 1.0 Release baseline (1bdbaa8c)
    "1a7edbb5"       # Optimise bessel_zeros
    "925fcb6d"       # Replace one powf with powi
    "06ef3897"       # Change to passing buffers, rather than allocating Vecs
    "de3dbb2c"       # Avoid re-calculation of consts in gamma_ln
    "7c4f2473"       # Tidy and optimise phase-factor calculation
    "47e6aa12"       # Tidy test_bessel_funcs
    "82482ad7"       # Improve naming of precision benchmarks in bessel-zeros
    "fc85169e"       # Optimise i_miller
    "200e40b6"       # Add new benchmarking crate
    "bb993e6f"       # Add dense and sequence workloads to benchmarks
    "0d148c90"       # Optimise i_power_series
    "2ad301c7"       # Optimise abs_ln_k calculations
    "208662dc"       # Remove string dependency by making bessel error details &'static str
    "49cb59e8"       # Clean up floating point calcs in right_half_plane
    "844fce94"       # Optimise recurrence
    "1b98b62f"       # Remove allocation in wronskian
    "aca0f1ab"       # Add fast path for integer order in Neumann normalisation loop
    "e0f60b22"       # Optimise horner polynomials and convolutions
    "51d7f3c6"       # Reduce tight-loop divisions
    "beca9b46"       # Hoist complex division in i_asymptotic
    "d3754470"       # Reduce calls to abs in i_ratios
    "b546b4ea"       # Improve register usage in i_power_series
    "85f2cece"       # Optimise large z miler seeds by adding fast path for real z
    "9ea47c70"       # Optimise loops in miller large seeds
    "HEAD"           # Current workspace
)

BENCH_COMMIT=$(git rev-parse perfomance-tweaks)
ORIGINAL_HEAD=$(git rev-parse HEAD)

cleanup() {
    echo "Restoring to original state: $ORIGINAL_HEAD"
    git checkout --quiet "$ORIGINAL_HEAD" 2>/dev/null || true
}
trap cleanup EXIT

for commit in "${COMMITS[@]}"; do
    echo "=================================================="
    echo " Benchmarking commit: $commit"
    echo "=================================================="

    # 1. Checkout the historical commit
    git checkout "$commit" --quiet

    # 2. Overlay the benchmark crate and workspace config from HEAD
    git checkout --quiet "$BENCH_COMMIT" -- bessel-benchmarks Cargo.toml Cargo.lock

    # 3. Run Criterion and save raw timing data under this commit name
    if cargo check -p bessel-benchmarks --bench workloads --features into --quiet 2>/dev/null; then
        cargo bench -p bessel-benchmarks --bench workloads --features into -- --save-baseline "$commit"
    else
        cargo bench -p bessel-benchmarks --bench workloads --no-default-features -- --save-baseline "$commit"
    fi

    # 4. Clean working tree before moving to the next commit
    git reset --hard --quiet
done

# Return to where you started
git checkout "$ORIGINAL_HEAD" --quiet
echo "All benchmarks completed!"
python3 bench_timeline.py
