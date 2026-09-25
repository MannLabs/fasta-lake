"""ESM-C representations of an explicitly recorded annotation-gap roster.

Model implementation: EvolutionaryScale ESM (now Biohub),
https://github.com/Biohub/esm. This adapter uses the pinned esm 3.2.3 ESMC
API and local published checkpoints, with FastaLake's explicit sequence
windowing and residue-weighted pooling. It does not train or redefine the model.
"""

from __future__ import annotations

import importlib.metadata
import json
import time
from pathlib import Path

from fasta_lake.multi_omics.exact import read_fasta, sha256, write_table


def sequence_windows(sequence, max_residues=1024):
    """Yield nonoverlapping full-coverage windows; never silently truncate a protein."""
    if type(max_residues) is not int or not 1 <= max_residues <= 2046:
        raise ValueError("max_residues must be an integer from 1 to 2046")
    if not sequence:
        raise ValueError("Cannot embed an empty protein sequence")
    for start in range(0, len(sequence), max_residues):
        yield start, sequence[start : start + max_residues]


def residue_embedding_sum(embeddings, length):
    """Sum residue embeddings in float64, excluding the ESM-C BOS/EOS tokens.

    A residue-weighted sum lets the caller combine a short final window with
    preceding windows without giving that shorter window disproportionate weight.
    """
    import numpy as np

    values = np.asarray(embeddings)
    if values.ndim != 2 or values.shape[0] != length + 2 or values.shape[1] < 1:
        raise ValueError("ESM-C output must contain one vector per residue plus BOS/EOS")
    if not np.isfinite(values).all():
        raise ValueError("ESM-C returned non-finite embeddings")
    return values[1:-1].sum(axis=0, dtype=np.float64)


def resolve_checkpoint(checkpoint, model):
    """Validate an explicit local checkpoint against its published upstream SHA256."""
    from fasta_lake.downloads import catalogue

    if model not in {"esmc_300m", "esmc_600m"}:
        raise ValueError("Choose esmc_300m or esmc_600m")
    resource = model.replace("_", "-")
    if checkpoint is None or not Path(checkpoint).is_file():
        raise FileNotFoundError(
            "Supply --checkpoint with a local ESM-C weight file. Download explicitly with "
            f"'fasta-lake resources download {resource} --directory /your/resources'."
        )
    path = Path(checkpoint).resolve(strict=True)
    expected = catalogue()[resource]["files"][0]["checksum"].removeprefix("sha256:")
    digest = sha256(path)
    if digest != expected:
        raise ValueError("ESM-C checkpoint differs from the published version for this model")
    return path, digest


def embed_unknown(
    annotation,
    output,
    *,
    model="esmc_300m",
    device="cpu",
    threads=2,
    max_residues=1024,
    checkpoint=None,
):
    """Embed every unknown from a completed, unchanged FastaLake annotation run.

    Unknown is defined by that run's recorded fields, not by a similarity score.
    ESM 3.2.3 with local 2024-12 ESM-C weights is the versioned adapter. Execution
    requires an explicitly supplied local checkpoint; no model download or remote
    inference API is used. CPU is the default; CUDA must be selected explicitly.

    Proteins longer than max_residues are covered by consecutive windows. Stored
    vectors average residue representations across all windows. This changes
    long-range context and must be matched when comparing reference embeddings.
    No functional annotation is inferred, and old prefix/all-token embeddings
    are not interchangeable with these vectors.
    """
    if model not in {"esmc_300m", "esmc_600m"}:
        raise ValueError("Choose esmc_300m or esmc_600m")
    if device not in {"cpu", "cuda"} or type(threads) is not int or threads < 1:
        raise ValueError("Choose cpu/cuda and a positive integer thread count")
    list(sequence_windows("A", max_residues))
    source, out = Path(annotation).resolve(), Path(output)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output exists: {out}; choose a new directory")
    receipt = json.loads((source / "COMPLETE.json").read_text())
    if receipt.get("status") != "PASS":
        raise ValueError("Embedding requires a completed annotation run")
    for name, expected in receipt["files"].items():
        if Path(name).name != name or sha256(source / name) != expected:
            raise ValueError("Annotation output changed before embedding")
    fasta = source / "unknown.fasta"
    fasta_hash = sha256(fasta)
    records = read_fasta(fasta) if fasta.stat().st_size else {}
    if len(records) != receipt["unknown_sequences"]:
        raise ValueError("Unknown FASTA and annotation receipt disagree")
    for query, (_, _, digest) in records.items():
        if query != "FLA_" + digest:
            raise ValueError("Unknown query ID does not match its exact sequence")
    if not records:
        out.mkdir(parents=True, exist_ok=False)
        result = {
            "status": "PASS",
            "embedded_sequences": 0,
            "backend_ran": False,
            "reason": "No unknowns under the recorded annotation rule",
            "unknown_fasta_sha256": fasta_hash,
        }
        (out / "COMPLETE.json").write_text(json.dumps(result, indent=2) + "\n")
        return result
    try:
        esm_version = importlib.metadata.version("esm")
    except importlib.metadata.PackageNotFoundError as error:
        raise RuntimeError(
            "ESM-C needs the optional embeddings environment: pip install '.[embeddings]'"
        ) from error
    if esm_version != "3.2.3":
        raise RuntimeError(f"This ESM-C adapter requires esm==3.2.3; found {esm_version}")
    import numpy as np
    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig
    from esm.tokenization import get_esmc_model_tokenizers

    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA was requested but is unavailable; choose --device cpu")
    torch.set_num_threads(threads)
    if torch.get_num_interop_threads() != 1:
        torch.set_num_interop_threads(1)
    model_size = "300" if model == "esmc_300m" else "600"
    checkpoint, checkpoint_hash = resolve_checkpoint(checkpoint, model)
    out.mkdir(parents=True, exist_ok=False)
    (out / "vectors").mkdir()
    provenance = {
        "esm": esm_version,
        "torch": torch.__version__,
        "model": model,
        "device": device,
        "checkpoint_sha256": checkpoint_hash,
        "checkpoint_path": str(checkpoint),
        "automatic_download": False,
        "threads": threads,
        "max_residues_per_window": max_residues,
        "window_overlap": 0,
        "pooling": "mean of all residue embeddings; BOS/EOS excluded; residue-weighted windows",
        "unknown_fields": receipt["unknown_fields"],
        "unknown_fasta_sha256": fasta_hash,
        "annotation_receipt_sha256": sha256(source / "COMPLETE.json"),
        "function_assignments_generated": False,
    }
    (out / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    started = time.monotonic()
    dimensions = (960, 15, 30) if model_size == "300" else (1152, 18, 36)
    with torch.device(device):
        client = ESMC(
            d_model=dimensions[0],
            n_heads=dimensions[1],
            n_layers=dimensions[2],
            tokenizer=get_esmc_model_tokenizers(),
            use_flash_attn=True,
        ).eval()
    state = torch.load(checkpoint, map_location=device, weights_only=True)
    client.load_state_dict(state)
    del state
    if device == "cuda":
        client = client.to(dtype=torch.bfloat16)
    manifest = []
    with torch.inference_mode():
        for query, (_, sequence, digest) in sorted(records.items()):
            total, windows = None, 0
            for _, window in sequence_windows(sequence, max_residues):
                encoded = client.encode(ESMProtein(sequence=window))
                logits = client.logits(
                    encoded, LogitsConfig(sequence=False, return_embeddings=True)
                )
                values = logits.embeddings.squeeze(0).float().cpu().numpy()
                subtotal = residue_embedding_sum(values, len(window))
                total = subtotal if total is None else total + subtotal
                windows += 1
            vector = (total / len(sequence)).astype(np.float32)
            expected_dimension = 960 if model_size == "300" else 1152
            if vector.shape != (expected_dimension,) or not np.isfinite(vector).all():
                raise ValueError("ESM-C output has invalid dimension or non-finite values")
            if np.linalg.norm(vector) == 0:
                raise ValueError("ESM-C returned a zero vector")
            relative = f"vectors/{digest}.npy"
            with (out / relative).open("xb") as stream:
                np.save(stream, vector, allow_pickle=False)
            manifest.append(
                {
                    "query": query,
                    "protein_sha256": digest,
                    "length": len(sequence),
                    "windows": windows,
                    "vector": relative,
                    "vector_sha256": sha256(out / relative),
                }
            )
            write_table(out / "vectors.tsv", list(manifest[0]), manifest)
    if sha256(fasta) != fasta_hash or sha256(checkpoint) != checkpoint_hash:
        raise ValueError("Unknown sequences or model checkpoint changed during embedding")
    result = {
        "status": "PASS",
        "embedded_sequences": len(manifest),
        "backend_ran": True,
        "seconds": time.monotonic() - started,
        "provenance": provenance,
        "vectors_tsv_sha256": sha256(out / "vectors.tsv"),
    }
    (out / "COMPLETE.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
