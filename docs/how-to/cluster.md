# Run on a cluster with Apptainer and Slurm

Most HPC systems do not allow Docker but do provide Apptainer (formerly
Singularity). Apptainer runs the same FastaLake image without root rights. This
page converts the image once, then runs a study as one Slurm job.

## 1. Convert the image once

Pick a published release from the
[releases page](https://github.com/MannLabs/fasta-lake/releases) and pull it into
a single `.sif` file on shared storage:

```bash
FASTALAKE_RELEASE=v1.1.0rc2          # an existing release tag
export APPTAINER_CACHEDIR=$PWD/.apptainer_cache
apptainer pull fasta-lake.sif "docker://ghcr.io/mannlabs/fasta-lake:$FASTALAKE_RELEASE"
apptainer exec --cleanenv fasta-lake.sif fasta-lake --version
```

Run `apptainer pull` on a node with internet access; the conversion needs a few
GB of temporary space, which `APPTAINER_CACHEDIR` keeps off your home quota.
A tested image downloaded from the Actions page can be converted offline with
`apptainer build fasta-lake.sif docker-archive://fastalake-docker-linux-amd64.tar.gz`.

!!! note "Images up to v1.1.0rc2"
    These images predate the `fasta-lake run` command. With them, replace
    `fasta-lake run` below by `python /opt/fastalake/tools/run_manifest.py`.
    The options are the same.

## 2. Bind your data and check the study

Apptainer sees only the directories you bind. Bind inputs read-only and give the
results a writable directory. `--cleanenv` keeps your login environment (Python
paths, conda variables) out of the container:

```bash
apptainer exec --cleanenv \
  --bind "$PWD/inputs:/data:ro" --bind "$PWD/results:/work" --pwd /work \
  fasta-lake.sif fasta-lake run \
  --manifest /data/samples.tsv --lake /data/lake.fasta \
  --lake-headers /data/joint_headers.tsv \
  --out /work/study_01 --sage /usr/local/bin/sage --threads 8 --validate-only
```

Paths inside the manifest are resolved inside the container, so write them
relative to the manifest (as in [the manifest guide](acquisition-manifest.md)) or
as `/data/...` paths. The Sage executable and the Rust tools are inside the image.

## 3. Submit the study as one job

The runner processes acquisitions one after another inside one job. Request as
many CPUs as `--threads` and size memory from a measured run of one complete
acquisition ([what computer do I need?](compute.md)); the memory request is
the hard limit Slurm enforces.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=fastalake_study_01
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.log
# Add --partition/--account as your site requires.
set -euo pipefail
apptainer exec --cleanenv \
  --bind "$PWD/inputs:/data:ro" --bind "$PWD/results:/work" --pwd /work \
  fasta-lake.sif fasta-lake run \
  --manifest /data/samples.tsv --lake /data/lake.fasta \
  --lake-headers /data/joint_headers.tsv \
  --out /work/study_01 --sage /usr/local/bin/sage \
  --threads "$SLURM_CPUS_PER_TASK"
```

Save it as `study_01.sbatch`, create `logs/` and `results/`, and submit with
`sbatch study_01.sbatch`. The output folder must not exist yet; a failed run
keeps its logs, so choose a new `--out` for the next attempt.

## 4. Check the result

A finished job is not yet a checked result. Confirm that `study_01/COMPLETE.json`
exists, that every stage in `commands.json` returned zero, and then open the files
listed in [read your results](read-outputs.md), starting with
`study_groups/annotated_study_group_matrix.tsv` and `study_groups/qc/QC_overview.png`.
