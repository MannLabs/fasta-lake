# Install FastaLake

For a first run, use Docker: it bundles the tools and example inputs. A cluster
that does not allow Docker can use the [source installation](#source-installation).

## Start with a working example

```bash
git clone https://github.com/MannLabs/fasta-lake.git
cd fasta-lake
bash tools/check_local_docker.sh docker_check_01
```

The first build downloads dependencies and compiles the tools, so allow several
minutes. The examples then run offline; a successful check ends with `PASS`.
Open [the demo guide](demo.md#open-these-first) to inspect your first results.

Prefer help from an agent or LLM? Use the [installation task file](../how-to/assistants.md#ready-to-copy-task-files).

## Choose another installation route

| Route | When to choose it | Needs |
|---|---|---|
| **Docker, built from the checkout** | Use Docker on a supported machine; examples run offline | Docker |
| **Tested image download** | You want the exact image a CI run tested, without building | Docker, Bash |
| **Source installation** | Development, or a cluster without Docker | Python 3.12, Rust toolchain, C compiler, Sage 0.14.6 |

The container includes Python, all nine FastaLake Rust executables, Sage 0.14.6, directLFQ 0.3.3, AlphaPeptTools 0.4.0 and the bundled public CAMPI acquisition windows. The complete Python runtime is pinned with distribution hashes in `requirements/container.txt` for Linux x86_64 and Python 3.12; the source route installs that same set.

## Docker, built from the checkout

From the repository root:

```bash
docker build --platform linux/amd64 -t fastalake:docker -f Dockerfile .
```

The first build downloads the base images and dependencies and compiles Rust. Allow several minutes and sufficient Docker memory for compilation. Builds default to two compiler jobs; use `--build-arg CARGO_BUILD_JOBS=1` if memory is limited. Subsequent builds reuse Docker's cache. Python wheels are mounted temporarily during installation so the runtime image does not retain a second copy of the dependency archives ([Docker build mounts](https://docs.docker.com/reference/dockerfile/#run---mounttypebind)).

The image targets **Linux x86_64**, matching the bundled Sage binary. On an Apple Silicon Mac, keep `--platform linux/amd64` in both build and run commands; Docker Desktop must support x86_64 emulation. The earlier core-only image was built and its synthetic, CAMPI and HSTOOL acceptance tests passed on an Apple Silicon Mac with Docker Desktop 4.46.0 (Engine 28.4.0) on 12 September 2026. The tests used amd64 emulation, a non-root user, no network and a read-only root filesystem. The exported archive was reloaded and its CLI verified.

## Download a tested image

Open a successful run in [Docker examples on GitHub](https://github.com/MannLabs/fasta-lake/actions/workflows/docker.yml)
and download its `fastalake-docker-linux-amd64` artifact. Choose the source commit
you want to use. GitHub keeps these artifacts for 14 days and applies repository
access rules. Unpack the download and open a terminal in that folder:

```bash
shasum -a 256 -c SHA256SUMS
docker load -i fastalake-docker-linux-amd64.tar.gz
docker tag fastalake:ci fastalake:docker
bash test_container.sh fastalake:docker demo_run_01
```

On Linux, `sha256sum --check SHA256SUMS` is an alternative. The download includes
the matching demo script, README, image identity and source commit, so this route
needs Docker and Bash without a source checkout or compiler. Keep those files
with the demo results. The workflow tests the image before export, then reloads
the archive and checks that the image identity matches.

## Versioned releases

The release workflow publishes a successfully tested image to
`ghcr.io/mannlabs/fasta-lake`, tagged with the release version and full source SHA.
Use the version shown on the corresponding GitHub release:

```bash
# Set this to an existing published release tag.
FASTALAKE_RELEASE=YOUR_RELEASE_TAG
docker pull --platform linux/amd64 "ghcr.io/mannlabs/fasta-lake:$FASTALAKE_RELEASE"
docker tag "ghcr.io/mannlabs/fasta-lake:$FASTALAKE_RELEASE" fastalake:docker
```

Until a versioned release and its container job complete, use a tested-image
download above or build from the checkout. Repository and package access
rules still apply; private packages require an authorized GHCR login. Published
non-prerelease versions also update `latest`; use a version or digest for a
reproducible run.

## Source installation

From the repository root, the locked route is Python 3.12 on Linux x86_64: it installs
the same hash-pinned package set as the container and CI
(`requirements/container.txt`), so a native run and a container run see the
same Python dependencies. One pinned package (directLFQ) builds from source, so
a C compiler is needed. Add the pinned Rust toolchain and a separately
installed Sage 0.14.6:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --require-hashes -r requirements/container.txt
python -m pip install --no-deps '.[analysis,directlfq]'
for crate in fasta_extractor_v2 parsimony_engine aggregation_engine; do
  cargo build --release --locked --manifest-path "rust/$crate/Cargo.toml"
done
export PATH="$PWD/rust/fasta_extractor_v2/target/release:$PWD/rust/parsimony_engine/target/release:$PWD/rust/aggregation_engine/target/release:$PATH"
```

The core package runs on Python 3.10 or newer; the analysis extras need 3.11 or newer; the tested and locked environment is 3.12.
Other Python versions and other platforms can install with
`python -m pip install '.[analysis,directlfq]'` instead of the two locked lines,
but then dependency versions float with PyPI and results are not the CI-tested
ones.

The Python wheel contains Python code; Rust executables are separate.
`FASTALAKE_BIN_DIR` can point to a directory containing those executables.
Run `python tools/run_manifest.py --help` and `fasta-lake --help` to see the
two entry points. The [function map](../reference/function-map.md) explains the library underneath.

## Verify the installation

The complete check builds the image and runs every bundled example inside it, offline, with a 4 GiB memory limit and two CPU slots:

```bash
bash tools/check_local_docker.sh docker_check_01
```

The command builds the combined Linux x86_64 image, then runs each test container
with a 4 GiB memory limit and two CPU slots, without network access and with a
read-only container filesystem. The build has separate resource requirements. Apple Silicon uses
Docker's x86_64 emulation. You do not need host Python, Rust, Sage, directLFQ
or AlphaPeptTools installations.

The first build downloads dependencies and compiles Rust. Tests cover the
synthetic selection example, measured public CAMPI windows, exact molecular
additions, real-data directLFQ/QC, and a separate synthetic PCA/missingness example.
The two-acquisition CAMPI example correctly skips PCA; no extra acquisitions
are invented for it.

Successful tests end with PASS. Keep the new `docker_check_01/` directory and
logs for the integration review. Existing output directories are refused.
The check builds locally and does not publish an image.

For actual inputs and registry release instructions, see
[the Docker guide](docker.md).

To verify an existing image, or carry it to a computer without internet:

```bash
bash tools/test_container.sh fastalake:docker docker_check
docker image inspect fastalake:docker --format '{{.Id}}'
docker save -o fastalake-docker-linux-amd64.tar fastalake:docker
gzip fastalake-docker-linux-amd64.tar
# On the receiving computer:
docker load -i fastalake-docker-linux-amd64.tar.gz
```

Run and verify the image before travel; then the example runs need no internet. GitHub Actions retains the tested image archive for 14 days on branch and pull-request runs. A published GitHub release additionally publishes that same tested image to GHCR. Record the source commit and image ID for reproducibility. Rust and all runtime Python dependencies are pinned. Base-image tags, OS packages and build tooling are not a fully immutable build lock. Record the image digest.

The container test runs all executable help checks, the synthetic examples, the public CAMPI real-data example, the additive metaG/metaT/metabolite workflow on real CAMPI spectra with synthetic molecular values, directLFQ and AlphaPeptTools QC on the CAMPI output, and PCA on a separately labelled three-acquisition synthetic fixture, with networking disabled and the container filesystem read-only. It also checks that Cargo and a C compiler are absent and that an existing output cannot be overwritten. Test results validate execution and expected outputs, not full-procedure FDR calibration.
