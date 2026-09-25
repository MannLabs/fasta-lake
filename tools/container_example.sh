#!/bin/sh
set -eu
case "${1:-}" in
    campi|hstool) cohort="$1"; shift ;;
    *) echo 'Usage: fastalake-example {campi|hstool} --out /work/NEW_DIRECTORY [--threads 2]' >&2
       exit 2 ;;
esac
if [[ ! -f "/opt/fastalake/examples/$cohort/inputs/samples.tsv" ]]; then
  echo "Example inputs are not distributed for $cohort; use the public campi example." >&2
  exit 2
fi
exec python "/opt/fastalake/examples/$cohort/run_example.py" "$@" \
    --runtime "${FASTALAKE_EXAMPLE_RUNTIME:-/opt/fastalake-runtime}"
