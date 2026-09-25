# syntax=docker/dockerfile:1@sha256:ecfaec9ed6d810b56388c508f4121597bfbba70d41a6dfeee4d8cad5f295fc32
# From the repository root:
# docker build --platform linux/amd64 -t fastalake:docker .
FROM rust:1.91.1-slim-bookworm@sha256:8514999d4786ef12efe89239e86b3d0a021b94b9d35108c8efe6c79ca7dc1a65 AS rust-builder
RUN test "$(dpkg --print-architecture)" = amd64 || \
    (echo 'The bundled Sage runtime requires --platform linux/amd64' >&2; exit 1)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential pkg-config libssl-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY rust/ ./rust/
COPY rust-toolchain.toml ./
ARG CARGO_BUILD_JOBS=2
ENV CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS}
RUN set -eux; \
    for crate in shared_peptide_classifier lake_stats parsimony_engine \
                 aggregation_engine fasta_export fasta_extractor_v2; do \
        cargo build --release --locked --manifest-path rust/$crate/Cargo.toml; \
    done; \
    mkdir -p /out/bin; \
    for crate in shared_peptide_classifier lake_stats parsimony_engine \
                 aggregation_engine fasta_export; do \
        cp rust/$crate/target/release/$crate /out/bin/; \
    done; \
    for bin in fasta_extractor lake_builder mass_kmer_anchor hash_fasta; do \
        cp rust/fasta_extractor_v2/target/release/$bin /out/bin/; \
    done

FROM python:3.12-slim-bookworm@sha256:782412e85d0f0984994c290652577d4018aff08145c85b262bb63dc0c7522254 AS python-builder
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY pyproject.toml README.md LICENSE MANIFEST.in ./
COPY requirements/ ./requirements/
COPY release/ ./release/
COPY fasta_lake/ ./fasta_lake/
RUN python -m pip wheel --no-cache-dir --wheel-dir /wheels \
    --require-hashes -r requirements/container.txt \
    && python -m pip wheel --no-cache-dir --no-deps --wheel-dir /wheels .

FROM python:3.12-slim-bookworm@sha256:782412e85d0f0984994c290652577d4018aff08145c85b262bb63dc0c7522254 AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates libssl3 libgomp1 \
    && rm -rf /var/lib/apt/lists/*
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba \
    FASTALAKE_BIN_DIR=/usr/local/bin \
    FASTALAKE_EXAMPLE_RUNTIME=/opt/fastalake-runtime
COPY --from=rust-builder /out/bin/ /usr/local/bin/
RUN --mount=type=bind,from=python-builder,source=/wheels,target=/wheels \
    python -m pip install --no-cache-dir --no-index --find-links=/wheels 'fasta-lake[directlfq,analysis]' \
    && python -m pip check
WORKDIR /opt/fastalake
# The examples verify these source files; targets and the toolchain are excluded.
COPY pyproject.toml rust-toolchain.toml README.md LICENSE ./
COPY requirements/ ./requirements/
COPY rust/ ./rust/
COPY fasta_lake/ ./fasta_lake/
COPY demo/ ./demo/
COPY examples/ ./examples/
COPY tools/run_manifest.py tools/group_study.py tools/prepare_container_runtime.py tools/check_analysis_example.py tools/check_annotation_contract.py ./tools/
COPY tools/container_example.sh /usr/local/bin/fastalake-example
RUN chmod 755 /usr/local/bin/fastalake-example \
    && python tools/prepare_container_runtime.py /opt/fastalake-runtime \
    && ln -s /opt/fastalake-runtime/bin/sage /usr/local/bin/sage \
    && chmod -R a+rX /opt/fastalake /opt/fastalake-runtime
RUN set -eux; \
    for bin in shared_peptide_classifier lake_stats parsimony_engine aggregation_engine \
               fasta_export fasta_extractor lake_builder mass_kmer_anchor hash_fasta sage; do \
        "$bin" --help >/dev/null; \
    done; \
    fasta-lake --help >/dev/null; \
    mkdir /work; chown 10001:10001 /work
USER 10001:10001
WORKDIR /work
CMD ["fasta-lake", "--help"]
