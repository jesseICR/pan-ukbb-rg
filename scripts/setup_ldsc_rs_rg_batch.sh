#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LDSC_RS_DIR="${LDSC_RS_DIR:-$ROOT_DIR/external/ldsc-rs-rg-batch}"
LDSC_RS_TARGET="${LDSC_RS_TARGET:-$ROOT_DIR/external/ldsc-rs-rg-batch-target}"
LDSC_RS_COMMIT="${LDSC_RS_COMMIT:-ddb8efb90f94b9ba4d15b13f1bf37a9f9fbdc68d}"
PATCH_FILE="$ROOT_DIR/patches/ldsc-rs-rg-batch.patch"
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-16}"

if [[ ! -f "$PATCH_FILE" ]]; then
  echo "error: missing patch file: $PATCH_FILE" >&2
  exit 1
fi

mkdir -p "$(dirname "$LDSC_RS_DIR")" "$(dirname "$LDSC_RS_TARGET")"

if [[ ! -d "$LDSC_RS_DIR/.git" ]]; then
  git clone https://github.com/sharifhsn/ldsc.git "$LDSC_RS_DIR"
fi

if git -C "$LDSC_RS_DIR" apply --reverse --check "$PATCH_FILE" >/dev/null 2>&1; then
  echo "ldsc-rs rg-batch patch already applied in $LDSC_RS_DIR"
else
  git -C "$LDSC_RS_DIR" fetch origin
  git -C "$LDSC_RS_DIR" checkout -q "$LDSC_RS_COMMIT"
  git -C "$LDSC_RS_DIR" reset --hard -q "$LDSC_RS_COMMIT"
  git -C "$LDSC_RS_DIR" apply "$PATCH_FILE"
fi

if command -v cargo >/dev/null 2>&1; then
  echo "Building patched ldsc-rs with local cargo..."
  CARGO_TARGET_DIR="$LDSC_RS_TARGET" cargo build --release --manifest-path "$LDSC_RS_DIR/Cargo.toml"
else
  if ! command -v docker >/dev/null 2>&1; then
    echo "error: neither cargo nor docker is available; install Rust or Docker to build patched ldsc-rs" >&2
    exit 1
  fi
  echo "Building patched ldsc-rs with Docker rust:1.91-bookworm..."
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -e CARGO_HOME=/work/external/.cargo-home \
    -e CARGO_TARGET_DIR=/work/external/ldsc-rs-rg-batch-target \
    -e CARGO_BUILD_JOBS="$CARGO_BUILD_JOBS" \
    -v "$ROOT_DIR:/work" \
    -w /work/external/ldsc-rs-rg-batch \
    rust:1.91-bookworm cargo build --release --jobs "$CARGO_BUILD_JOBS"
fi

"$LDSC_RS_TARGET/release/ldsc" --help >/dev/null
echo "Built patched ldsc-rs: $LDSC_RS_TARGET/release/ldsc"
