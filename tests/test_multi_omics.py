"""Tests for the FL_MO multi-omics evidence gate (Methods §M5.5).

Covers the gate truth-table (OR/AND), the per-sample top-1% metaG threshold
with tie-to-smallest handling, the file-level Stage-3 runner round-trip, and
the reproduction of the published razor per-sample median (13,776).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta_lake.hashing import sequence_hash
from fasta_lake.multi_omics import (
    Gate,
    build_evidence_lookup,
    gate_members,
    passes_tpm,
    run_fl_mo_sample,
    top1_threshold,
)

# --------------------------------------------------------------------------
# top-1% threshold
# --------------------------------------------------------------------------


def test_top1_threshold_basic():
    # 100 non-zero values 1..100; top 1% is the single largest value (100).
    vals = list(range(1, 101))
    assert top1_threshold(vals, 1.0) == 100


def test_top1_threshold_ignores_zeros():
    vals = [0.0, 0.0, 0.0] + list(range(1, 101))
    # zeros are excluded from the distribution, so cutoff unchanged
    assert top1_threshold(vals, 1.0) == 100


def test_top1_threshold_ties_include_all():
    # boundary value tied: cutoff must drop to the smallest tied value so all
    # tied proteins are admitted (Methods §M5.5).
    vals = [1, 2, 3, 4, 5, 5, 5, 5, 5, 5]  # 10 values; top ~1% lands on a 5
    cut = top1_threshold(vals, 10.0)  # top 10% = 1 entry, but boundary tied at 5
    assert cut == 5
    # every tied 5 passes
    assert all(v >= cut for v in [5, 5, 5])


def test_top1_threshold_empty():
    assert top1_threshold([], 1.0) == 0.0
    assert top1_threshold([0.0, 0.0], 1.0) == 0.0


def test_passes_tpm():
    # metaG passes when >= cutoff and cutoff > 0
    assert passes_tpm(10.0, 0.0, metag_cutoff=5.0) == (True, False)
    assert passes_tpm(4.0, 0.0, metag_cutoff=5.0) == (False, False)
    # metaT passes strictly above 1.0
    assert passes_tpm(0.0, 1.5, metag_cutoff=5.0) == (False, True)
    assert passes_tpm(0.0, 1.0, metag_cutoff=5.0) == (False, False)
    # cutoff 0 => no metaG rescue possible
    assert passes_tpm(100.0, 0.0, metag_cutoff=0.0) == (False, False)


# --------------------------------------------------------------------------
# gate truth-table (the ARCHITECTURE.md fixture cases)
# --------------------------------------------------------------------------


def _gate_setup():
    """Construct inputs exercising each evidence combination.

    Cluster reps:
      REP_DN     - de novo only (no TPM)        -> OR yes, AND no
      REP_MG     - metaG-top1 only (no de novo) -> OR yes, AND no
      REP_MT     - metaT>1 only (no de novo)    -> OR yes, AND no
      REP_BOTH   - de novo AND metaG-top1       -> OR yes, AND yes
      REP_LOWTPM - de novo, low TPM             -> OR yes (denovo), AND no
    """
    denovo = {"REP_DN", "REP_BOTH", "REP_LOWTPM"}

    # assembly proteins (id -> sequence). Distinct sequences -> distinct hashes.
    seqs = {
        "k141_1": "MAAAAAAAAAK",  # -> REP_MG   (high metaG)
        "k141_2": "MCCCCCCCCCK",  # -> REP_MT   (high metaT)
        "k141_3": "MDDDDDDDDDK",  # -> REP_BOTH (high metaG)
        "k141_4": "MEEEEEEEEEK",  # -> REP_LOWTPM (low TPM, won't pass)
        "k141_5": "MFFFFFFFFFK",  # high metaG -> REP_DN (de novo, also TPM)?
    }
    hash_to_rep = {
        sequence_hash(seqs["k141_1"]): "REP_MG",
        sequence_hash(seqs["k141_2"]): "REP_MT",
        sequence_hash(seqs["k141_3"]): "REP_BOTH",
        sequence_hash(seqs["k141_4"]): "REP_LOWTPM",
        sequence_hash(seqs["k141_5"]): "REP_DN",
    }
    # metaG distribution: many low values + a few high so top1% cutoff is high.
    # Build 100 low entries plus the high signal entries.
    tpm_rows = [(f"bg_{i}", 1.0, 0.0) for i in range(100)]
    tpm_rows += [
        ("k141_1", 999.0, 0.0),  # metaG high -> passes metaG
        ("k141_2", 0.0, 5.0),  # metaT 5 > 1 -> passes metaT
        ("k141_3", 999.0, 0.0),  # metaG high -> passes metaG (REP_BOTH)
        ("k141_4", 1.0, 0.5),  # low metaG, metaT<=1 -> fails
        # k141_5 intentionally absent from TPM (REP_DN has no TPM evidence)
    ]
    return denovo, tpm_rows, seqs, hash_to_rep


def test_or_and_gate_truth_table():
    denovo, tpm_rows, seqs, htr = _gate_setup()
    rows, stats = build_evidence_lookup(denovo, tpm_rows, seqs, htr)
    by_id = {r.cluster_rep_id: r for r in rows}

    # de novo only
    assert by_id["REP_DN"].has_denovo and not by_id["REP_DN"].has_tpm
    # metaG-rescue only
    assert by_id["REP_MG"].source == "TPM_RESCUE"
    assert by_id["REP_MG"].has_metaG_top1 and not by_id["REP_MG"].has_denovo
    # metaT-rescue only
    assert by_id["REP_MT"].has_metaT_1 and not by_id["REP_MT"].has_denovo
    # de novo AND metaG
    assert by_id["REP_BOTH"].has_denovo and by_id["REP_BOTH"].has_metaG_top1
    # low TPM de novo: included as DENOVO but no TPM flags
    assert by_id["REP_LOWTPM"].has_denovo and not by_id["REP_LOWTPM"].has_tpm

    or_set = gate_members(rows, Gate.OR)
    and_set = gate_members(rows, Gate.AND)

    # OR-gate admits everything with any evidence
    assert or_set == {"REP_DN", "REP_MG", "REP_MT", "REP_BOTH", "REP_LOWTPM"}
    # AND-gate admits only reps with de novo AND TPM evidence
    assert and_set == {"REP_BOTH"}

    # stats sanity
    assert stats["or_gate_members"] == 5
    assert stats["and_gate_members"] == 1
    assert stats["denovo_reps"] == 3


def test_low_tpm_low_denovo_excluded():
    # a rep with neither de novo nor passing TPM never appears / never admitted
    denovo = {"REP_DN"}
    seqs = {"k141_x": "MAAAAAAAAAK"}
    htr = {sequence_hash(seqs["k141_x"]): "REP_X"}
    tpm_rows = [("k141_x", 0.5, 0.5)] + [(f"bg_{i}", 1.0, 0.0) for i in range(100)]
    rows, _ = build_evidence_lookup(denovo, tpm_rows, seqs, htr)
    ids = {r.cluster_rep_id for r in rows}
    assert "REP_X" not in ids
    assert gate_members(rows, Gate.OR) == {"REP_DN"}


def test_threshold_sensitivity():
    # widening the metaG top-pct admits more metaG-rescued reps
    denovo = set()
    seqs = {f"k141_{i}": f"M{'A' * i}K" for i in range(2, 12)}
    htr = {sequence_hash(s): f"REP_{i}" for i, s in zip(range(2, 12), seqs.values())}
    # 100 background lows + the 10 candidates with ascending metaG
    tpm_rows = [(f"bg_{i}", 0.1, 0.0) for i in range(100)]
    tpm_rows += [(f"k141_{i}", float(i * 100), 0.0) for i in range(2, 12)]
    rows_strict, _ = build_evidence_lookup(denovo, tpm_rows, seqs, htr, metag_top_pct=1.0)
    rows_loose, _ = build_evidence_lookup(denovo, tpm_rows, seqs, htr, metag_top_pct=10.0)
    assert len(gate_members(rows_loose, Gate.OR)) >= len(gate_members(rows_strict, Gate.OR))


# --------------------------------------------------------------------------
# file-level runner round-trip
# --------------------------------------------------------------------------


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def test_runner_roundtrip_or_gate(tmp_path):
    # de-novo FASTA with two reps (multi-line sequence for one)
    denovo_fa = tmp_path / "denovo.fasta"
    _write(denovo_fa, ">REP_DN desc here\nMKKKKKKKKK\nMKKKKK\n>REP_BOTH\nMDDDDDDDDDK\n")

    # assembly FASTA (MEGAHIT predicted genes)
    asm_fa = tmp_path / "asm.fasta"
    asm_seqs = {"k141_1": "MAAAAAAAAAK", "k141_3": "MDDDDDDDDDK"}
    _write(asm_fa, "".join(f">{k}\n{v}\n" for k, v in asm_seqs.items()))

    # cluster rep FASTA (source of verbatim rescued sequences)
    rep_fa = tmp_path / "reps.fasta"
    _write(
        rep_fa,
        ">REP_DN x\nMKKKKKKKKKMKKKKK\n>REP_BOTH\nMDDDDDDDDDK\n"
        ">REP_MG some annotation\nMAAAAAAAAAK\n",
    )

    # hash_to_rep
    htr = tmp_path / "h2r.tsv"
    _write(
        htr,
        "hash\trep\n"
        f"{sequence_hash(asm_seqs['k141_1'])}\tREP_MG\n"
        f"{sequence_hash(asm_seqs['k141_3'])}\tREP_BOTH\n",
    )

    # TPM: background lows + k141_1 high metaG (-> REP_MG rescue), k141_3 high
    tpm = tmp_path / "tpm.tsv"
    lines = ["prodigal_protein\tmetaG_tpm\tmetaT_tpm"]
    lines += [f"bg_{i}\t1.0\t0.0" for i in range(100)]
    lines += ["k141_1\t999.0\t0.0", "k141_3\t999.0\t0.0"]
    _write(tpm, "\n".join(lines) + "\n")

    out_fa = tmp_path / "out.fasta"
    out_lk = tmp_path / "out_lookup.tsv"
    out_st = tmp_path / "out.stats.json"

    stats = run_fl_mo_sample(
        sample="S1",
        denovo_fasta=str(denovo_fa),
        tpm_tsv=str(tpm),
        assembly_fasta=str(asm_fa),
        hash_to_rep=str(htr),
        cluster_rep_fasta=str(rep_fa),
        output_fasta=str(out_fa),
        output_lookup=str(out_lk),
        output_stats=str(out_st),
        gate=Gate.OR,
    )

    # OR-gate: both de novo reps + the metaG-rescued REP_MG
    headers = [
        line[1:].split()[0] for line in out_fa.read_text().splitlines() if line.startswith(">")
    ]
    assert set(headers) == {"REP_DN", "REP_BOTH", "REP_MG"}
    assert stats["fasta_denovo_written"] == 2
    assert stats["fasta_rescued_written"] == 1
    # rescued header is copied verbatim from rep FASTA (tagged), seq preserved
    assert "MAAAAAAAAAK" in out_fa.read_text()


def test_runner_and_gate_subset_of_or(tmp_path):
    denovo_fa = tmp_path / "denovo.fasta"
    _write(denovo_fa, ">REP_DN\nMKKKKKKKKK\n>REP_BOTH\nMDDDDDDDDDK\n")
    asm_fa = tmp_path / "asm.fasta"
    _write(asm_fa, ">k141_3\nMDDDDDDDDDK\n>k141_1\nMAAAAAAAAAK\n")
    rep_fa = tmp_path / "reps.fasta"
    _write(rep_fa, ">REP_DN\nMKKKKKKKKK\n>REP_BOTH\nMDDDDDDDDDK\n>REP_MG\nMAAAAAAAAAK\n")
    htr = tmp_path / "h2r.tsv"
    _write(
        htr, f"{sequence_hash('MDDDDDDDDDK')}\tREP_BOTH\n{sequence_hash('MAAAAAAAAAK')}\tREP_MG\n"
    )
    tpm = tmp_path / "tpm.tsv"
    lines = ["protein\tmetaG_tpm\tmetaT_tpm"]
    lines += [f"bg_{i}\t1.0\t0.0" for i in range(100)]
    lines += ["k141_3\t999.0\t0.0", "k141_1\t999.0\t0.0"]
    _write(tpm, "\n".join(lines) + "\n")

    out_fa = tmp_path / "and.fasta"
    stats = run_fl_mo_sample(
        sample="S1",
        denovo_fasta=str(denovo_fa),
        tpm_tsv=str(tpm),
        assembly_fasta=str(asm_fa),
        hash_to_rep=str(htr),
        cluster_rep_fasta=str(rep_fa),
        output_fasta=str(out_fa),
        output_lookup=str(tmp_path / "lk.tsv"),
        gate=Gate.AND,
    )
    headers = {
        line[1:].split()[0] for line in out_fa.read_text().splitlines() if line.startswith(">")
    }
    # AND-gate keeps only REP_BOTH (de novo AND metaG); REP_DN dropped (no TPM)
    assert headers == {"REP_BOTH"}
    assert stats["and_gate_members"] == 1


# --------------------------------------------------------------------------
# reproduce the published 13,776
# --------------------------------------------------------------------------


def test_reproduce_published_median():
    from pathlib import Path

    from fasta_lake.multi_omics import reproduce as repro_mod

    # only run when the committed verification TSV is present on this checkout
    if not Path(repro_mod.FL_MO_RAZOR_TSV).exists():
        import pytest

        pytest.skip("verification TSV not present in this checkout")
    r = repro_mod.reproduce(verbose=False)
    assert r["fl_mo_razor_median"] == 13776
    assert r["matches_published"]


# --------------------------------------------------------------------------
# byte-identical outputs across PYTHONHASHSEED (set-order regression)
# --------------------------------------------------------------------------

_HASHSEED_SCRIPT = r"""
import sys
from fasta_lake.hashing import sequence_hash
from fasta_lake.multi_omics import Gate, run_fl_mo_sample
out = sys.argv[1]
n = 12
denovo = "".join(f">REP_DN{i:02d} d\nM{'K' * (i + 5)}\n" for i in range(n))
asm = {f"k141_{i}": f"M{'A' * (i + 5)}K" for i in range(n)}
open(f"{out}/denovo.fasta", "w").write(denovo)
open(f"{out}/asm.fasta", "w").write("".join(f">{k}\n{v}\n" for k, v in asm.items()))
open(f"{out}/reps.fasta", "w").write(
    denovo + "".join(f">REP_MG{i:02d} ann\n{asm[f'k141_{i}']}\n" for i in range(n))
)
open(f"{out}/h2r.tsv", "w").write(
    "hash\trep\n"
    + "".join(f"{sequence_hash(v)}\tREP_MG{i:02d}\n" for i, v in enumerate(asm.values()))
)
tpm = ["prodigal_protein\tmetaG_tpm\tmetaT_tpm"] + [f"bg_{i}\t1.0\t0.0" for i in range(100)]
tpm += [f"k141_{i}\t999.0\t0.0" for i in range(n)]
open(f"{out}/tpm.tsv", "w").write("\n".join(tpm) + "\n")
run_fl_mo_sample(
    sample="S", denovo_fasta=f"{out}/denovo.fasta", tpm_tsv=f"{out}/tpm.tsv",
    assembly_fasta=f"{out}/asm.fasta", hash_to_rep=f"{out}/h2r.tsv",
    cluster_rep_fasta=f"{out}/reps.fasta", output_fasta=f"{out}/out.fasta",
    output_lookup=f"{out}/out_lookup.tsv", output_stats=None, gate=Gate.OR,
)
"""


def test_outputs_byte_identical_across_hash_seeds(tmp_path):
    """evidence_lookup.tsv and the rescued FASTA must not depend on set iteration order."""
    import subprocess

    script = tmp_path / "run.py"
    script.write_text(_HASHSEED_SCRIPT)
    outputs = []
    for seed in ("0", "1", "4242"):
        out = tmp_path / f"seed{seed}"
        out.mkdir()
        env = dict(os.environ, PYTHONHASHSEED=seed)
        subprocess.run([sys.executable, str(script), str(out)], check=True, env=env)
        outputs.append(((out / "out_lookup.tsv").read_bytes(), (out / "out.fasta").read_bytes()))
    assert outputs[0] == outputs[1] == outputs[2]
    lookup_ids = [ln.split("\t")[0] for ln in outputs[0][0].decode().splitlines()[1:]]
    assert lookup_ids == sorted(lookup_ids[:12]) + sorted(lookup_ids[12:])
