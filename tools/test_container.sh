#!/usr/bin/env bash
# Test an already-built image. Keep the output directory as the validation record.
set -euo pipefail
image="${1:-fastalake:docker}"
out="${2:-docker_check}"
if [[ -e "$out" ]]; then
  echo "Output exists: $out; choose a new directory" >&2
  exit 2
fi
mkdir -p "$out"
out="$(cd "$out" && pwd)"
run() {
  docker run --rm --platform linux/amd64 --network none --read-only \
    --memory 4g --memory-swap 4g --cpus 2 \
    --tmpfs /tmp:rw,nosuid,nodev --user "$(id -u):$(id -g)" \
    --mount "type=bind,source=$out,target=/work" "$image" "$@"
}
docker image inspect "$image" > "$out/image.json"
run sh -ec '
  for compiler in cargo rustc gcc; do
    if command -v "$compiler" >/dev/null 2>&1; then
      echo "Unexpected build tool in runtime: $compiler" >&2
      exit 1
    fi
  done
  for bin in shared_peptide_classifier lake_stats parsimony_engine aggregation_engine \
             fasta_export fasta_extractor lake_builder mass_kmer_anchor hash_fasta sage fasta-lake; do
    "$bin" --help >/dev/null
  done
' > "$out/runtime.log" 2>&1
run bash /opt/fastalake/demo/run_demo.sh /work/demo > "$out/demo.log" 2>&1
run python -c '
from pathlib import Path
p = Path("/work/demo")
count = lambda f: sum(line.startswith(">") for line in f.open())
assert count(p / "evidence_lake.fasta") == 14
for sample, expected in zip(("SAMPLE_A", "SAMPLE_B", "SAMPLE_C"), (4, 4, 3)):
    assert count(p / sample / "inference" / (sample + "_razor.fasta")) == expected
' >> "$out/demo.log" 2>&1
for cohort in campi; do
  run fastalake-example "$cohort" --out "/work/$cohort" --threads 2 > "$out/$cohort.log" 2>&1
done
run python /opt/fastalake/examples/molecular/run_example.py --out /work/molecular > "$out/molecular.log" 2>&1
run fasta-lake quantify-study --groups /work/campi/study --method directlfq \
  --out /work/campi_directlfq --threads 1 > "$out/directlfq.log" 2>&1
run fasta-lake analyze-study --groups /work/campi/study \
  --quantification /work/campi_directlfq --out /work/campi_analysis > "$out/analysis.log" 2>&1
run python /opt/fastalake/tools/check_analysis_example.py \
  /work/campi/study /work/campi_directlfq /work/campi_analysis >> "$out/analysis.log" 2>&1
run python /opt/fastalake/examples/analysis/run_example.py --out /work/synthetic_analysis > "$out/synthetic_analysis.log" 2>&1
run python /opt/fastalake/examples/multiomics/run_example.py --campi /work/campi \
  --out /work/multiomics --threads 2 --reference-chunk-mib 1 --laptop > "$out/multiomics.log" 2>&1
run python /opt/fastalake/examples/chunked/run_example.py --campi /work/campi \
  --out /work/chunked --laptop > "$out/chunked.log" 2>&1
# Keep all six search modes in the 4 GiB offline teaching check. The separate
# host acceptance also tests 7-30 aa and measures the larger unspecific index.
run python /opt/fastalake/examples/searches/run_example.py --campi /work/campi \
  --out /work/search_choices --threads 2 --search-max-length 15 > "$out/search_choices.log" 2>&1
run python /opt/fastalake/examples/functional/run_example.py --out /work/functional > "$out/functional.log" 2>&1
run python /opt/fastalake/tools/check_annotation_contract.py --campi /work/campi \
  --out /work/annotation_contract > "$out/annotation_contract.log" 2>&1
run python /opt/fastalake/examples/community/run_example.py --groups /work/campi/study \
  --out /work/community > "$out/community.log" 2>&1
run fasta-lake resources list > "$out/resource_catalogue.log" 2>&1
run fasta-lake annotate --help > "$out/annotation_help.log" 2>&1
if run fastalake-example campi --out /work/campi --threads 2 > "$out/refuse_existing.log" 2>&1; then
  echo 'An existing output was unexpectedly accepted' >&2
  exit 1
fi
grep -q 'Output already exists' "$out/refuse_existing.log"
printf 'PASS: prebuilt runtime, synthetic demo, CAMPI, molecular input, directLFQ, AlphaPeptTools, additive multi-omics, search options, experimental functional profiles, Unipept taxonomy, evidence networks, eggNOG joins and output protection.\n'
