"""Resource planning must respect parent limits and retain laptop headroom."""

import importlib.util
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from fasta_lake import capacity
from fasta_lake.cli import cli

GIB = 1024**3


def machine(total=32, available=24, cpus=8):
    return {
        "memory_capacity_bytes": total * GIB,
        "memory_available_bytes": available * GIB,
        "cpu_slots": cpus,
    }


def test_32_gib_machine_uses_shared_threads_and_leaves_headroom():
    plan = capacity.plan_laptop(capacity=machine(), memory_gib=11)
    assert plan["memory_budget_bytes"] == 11 * GIB
    assert plan["threads_per_process"] == 4
    assert plan["parallel_processes"] == 1
    assert plan["reference_chunk_mib"] == 256
    assert plan["blas_threads"] == 1


def test_busy_machine_and_fractional_cpu_still_have_a_valid_small_plan():
    plan = capacity.plan_laptop(capacity=machine(32, 2, 1))
    assert plan["memory_budget_bytes"] == int(1.6 * GIB)
    assert plan["threads_per_process"] == 1
    assert plan["reference_chunk_mib"] == 102


@pytest.mark.parametrize(
    "kwargs",
    [
        {"memory_gib": 30},
        {"memory_gib": float("nan")},
        {"memory_gib": 0},
        {"memory_fraction": 0},
        {"memory_fraction": 1.1},
        {"cpu_fraction": float("inf")},
        {"threads": 9},
        {"threads": 0},
        {"threads": True},
        {"chunk_mib": 1000},
    ],
)
def test_invalid_or_excessive_budgets_are_rejected(kwargs):
    with pytest.raises(ValueError):
        capacity.plan_laptop(capacity=machine(), **kwargs)


def test_no_memory_headroom_fails():
    with pytest.raises(ValueError, match="Insufficient"):
        capacity.plan_laptop(capacity=machine(32, 0, 8))


def cgroup_fixture(tmp_path, *, version=2):
    proc = tmp_path / "proc"
    (proc / "self").mkdir(parents=True)
    mount = tmp_path / "cgroup"
    leaf = mount / "job/step/task"
    leaf.mkdir(parents=True)
    (proc / "meminfo").write_text("MemTotal: 268435456 kB\nMemAvailable: 209715200 kB\n")
    if version == 2:
        (proc / "self/cgroup").write_text("0::/job/step/task\n")
        (proc / "self/mountinfo").write_text(f"38 36 0:25 / {mount} rw - cgroup2 cgroup rw\n")
        (leaf / "memory.max").write_text("max\n")
        (mount / "job/memory.max").write_text(str(12 * GIB))
        (mount / "job/memory.current").write_text(str(8 * GIB))
        (mount / "job/memory.stat").write_text(f"inactive_file {2 * GIB}\n")
        (mount / "job/cpu.max").write_text("250000 100000\n")
    else:
        (proc / "self/cgroup").write_text("4:memory,cpu:/job/step/task\n")
        (proc / "self/mountinfo").write_text(
            f"38 36 0:25 / {mount} rw - cgroup cgroup rw,memory,cpu\n"
        )
        (leaf / "memory.limit_in_bytes").write_text(str(2**63 - 4096))
        (mount / "job/memory.limit_in_bytes").write_text(str(12 * GIB))
        (mount / "job/memory.usage_in_bytes").write_text(str(8 * GIB))
        (mount / "job/memory.stat").write_text(f"total_inactive_file {2 * GIB}\n")
        (mount / "job/cpu.cfs_quota_us").write_text("250000")
        (mount / "job/cpu.cfs_period_us").write_text("100000")
    return proc, mount


@pytest.mark.parametrize("version", [1, 2])
def test_limits_include_cgroup_ancestors_and_reclaimable_cache(tmp_path, monkeypatch, version):
    proc, _ = cgroup_fixture(tmp_path, version=version)
    monkeypatch.setattr(capacity.sys, "platform", "linux")
    monkeypatch.setattr(capacity.os, "cpu_count", lambda: 128)
    monkeypatch.setattr(capacity.os, "sched_getaffinity", lambda _: {0, 1, 2, 3}, raising=False)
    detected = capacity.detect_capacity(
        proc=proc, environ={"SLURM_MEM_PER_NODE": "8192", "SLURM_CPUS_PER_TASK": "4"}
    )
    assert detected["memory_capacity_bytes"] == 8 * GIB
    assert detected["memory_available_bytes"] == 6 * GIB
    assert detected["cpu_slots"] == 2
    assert len(detected["cgroup_observations"]) == 1


def test_private_cgroup_namespace_reads_its_mounted_root(tmp_path):
    proc, mount = cgroup_fixture(tmp_path)
    (proc / "self/cgroup").write_text("0::/\n")
    (proc / "self/mountinfo").write_text(
        f"38 36 0:25 /docker/container {mount} rw - cgroup2 cgroup rw\n"
    )
    (mount / "memory.max").write_text(str(4 * GIB))
    assert capacity._cgroup_limits(proc)[0]["memory_limit_bytes"] == 4 * GIB


def test_slurm_per_cpu_and_affinity_limits(tmp_path, monkeypatch):
    proc, _ = cgroup_fixture(tmp_path)
    (proc / "self/mountinfo").write_text("")
    monkeypatch.setattr(capacity.sys, "platform", "linux")
    monkeypatch.setattr(capacity.os, "cpu_count", lambda: 128)
    monkeypatch.setattr(capacity.os, "sched_getaffinity", lambda _: {3}, raising=False)
    detected = capacity.detect_capacity(
        proc=proc, environ={"SLURM_MEM_PER_CPU": "1024", "SLURM_CPUS_PER_TASK": "2"}
    )
    assert detected["memory_capacity_bytes"] == 2 * GIB
    assert detected["cpu_slots"] == 1


def test_macos_memory_reads_page_size_and_available_pages(monkeypatch):
    monkeypatch.setattr(capacity.sys, "platform", "darwin")
    monkeypatch.setattr(
        capacity.subprocess,
        "check_output",
        lambda argv, **kw: (
            str(32 * GIB)
            if argv[0] == "sysctl"
            else (
                "Mach Virtual Memory Statistics: (page size of 16384 bytes)\n"
                "Pages free: 100.\nPages inactive: 200.\nPages speculative: 50.\n"
            )
        ),
    )
    assert capacity._host_memory(Path("/missing")) == (32 * GIB, 350 * 16384)


def test_installed_command_reports_plan_without_creating_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(capacity, "detect_capacity", machine)
    result = CliRunner().invoke(cli, ["plan-resources", "--memory-budget-gib", "11"])
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["memory_budget_bytes"] == 11 * GIB
    assert list(tmp_path.iterdir()) == []


def test_manifest_preserves_explicit_zero_thread_rejection(tmp_path):
    path = Path(__file__).resolve().parents[1] / "tools/run_manifest.py"
    spec = importlib.util.spec_from_file_location("capacity_runner", path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    with pytest.raises(SystemExit) as error:
        runner.main(
            [
                "--manifest",
                str(tmp_path / "missing"),
                "--out",
                str(tmp_path / "out"),
                "--threads",
                "0",
            ]
        )
    assert error.value.code == 2
    assert not (tmp_path / "out").exists()


def test_laptop_helper_limits_are_set_before_preflight_imports(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[1] / "tools/run_manifest.py"
    spec = importlib.util.spec_from_file_location("laptop_preflight_runner", path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    monkeypatch.setattr(capacity, "detect_capacity", machine)
    names = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")
    for name in names:
        monkeypatch.setenv(name, "8")

    def preflight(_):
        assert all(capacity.os.environ[name] == "1" for name in names)
        raise ValueError("reached preflight after applying helper limits")

    monkeypatch.setattr(runner, "load_manifest", preflight)
    lake = tmp_path / "lake.fasta"
    lake.write_text(">test\nMPEPTIDEK\n")
    with pytest.raises(ValueError, match="reached preflight"):
        runner.main(
            [
                "--manifest",
                str(tmp_path / "manifest.tsv"),
                "--lake",
                str(lake),
                "--out",
                str(tmp_path / "out"),
                "--laptop",
            ]
        )
    assert not (tmp_path / "out").exists()
