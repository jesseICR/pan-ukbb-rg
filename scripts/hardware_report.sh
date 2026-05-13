#!/usr/bin/env bash
set -euo pipefail

echo "# Hardware"
date --iso-8601=seconds
echo
echo "## CPU"
lscpu
echo
echo "## Memory"
free -h
echo
echo "## Disk"
df -h .
echo
echo "## GPU"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true

