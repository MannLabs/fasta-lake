#!/usr/bin/env bash
# From a clone or extracted source bundle; requires only Docker and Bash.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
out="${1:-docker_check_01}"
image="fastalake:local-check"
if [[ -e "$out" || -L "$out" ]]; then
  echo "Output exists: $out; choose a new directory" >&2
  exit 2
fi
docker version >/dev/null
docker build --platform linux/amd64 -t "$image" -f "$repo/Dockerfile" "$repo"
exec bash "$repo/tools/test_container.sh" "$image" "$out"
