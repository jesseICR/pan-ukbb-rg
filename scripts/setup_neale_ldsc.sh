#!/usr/bin/env bash
set -euo pipefail

repo_url="https://github.com/astheeggeggs/ldsc.git"
commit="a4ee4c8aa065a1c9a586c3b678e9b3040bbebafc"
dest="${LDSC_DIR:-external/ldsc-neale}"

mkdir -p external

if [[ ! -d "$dest/.git" ]]; then
  git clone "$repo_url" "$dest"
fi

git -C "$dest" fetch --all --tags
git -C "$dest" checkout "$commit"

echo "Neale LDSC ready at $dest ($commit)"
