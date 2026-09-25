"""Conservative laptop planning within visible host, container and scheduler limits."""

from __future__ import annotations

import math
import os
import re
import subprocess
import sys
from pathlib import Path

GIB = 1024**3
MIB = 1024**2


def _integer(text):
    """Parse a positive capacity below 2**60, or return None for an unusable limit."""
    try:
        value = int(text)
        return value if 0 < value < 2**60 else None
    except (ValueError, TypeError):
        return None


def _read(path):
    """Read a stripped system-limit file, returning an empty string on I/O failure."""
    try:
        return path.read_text().strip()
    except OSError:
        return ""


def _host_memory(proc):
    """Return total and available host memory in bytes from macOS or procfs."""
    if sys.platform == "darwin":
        total = int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True))
        vm = subprocess.check_output(["vm_stat"], text=True)
        page_size = int(re.search(r"page size of (\d+) bytes", vm)[1])
        fields = dict(re.findall(r"^([^:\n]+):\s+(\d+)\.", vm, re.MULTILINE))
        available = (
            sum(
                int(fields.get(n, 0))
                for n in (
                    "Pages free",
                    "Pages inactive",
                    "Pages speculative",
                )
            )
            * page_size
        )
        return total, min(total, available)
    fields = {}
    for line in _read(proc / "meminfo").splitlines():
        key, _, value = line.partition(":")
        fields[key] = int(value.split()[0]) * 1024
    if "MemTotal" not in fields or "MemAvailable" not in fields:
        raise ValueError("Cannot detect available RAM; laptop planning supports Linux and macOS")
    return fields["MemTotal"], fields["MemAvailable"]


def _unescape_mount(path):
    """Decode octal escapes in a Linux mount path."""
    return re.sub(r"\\([0-7]{3})", lambda m: chr(int(m[1], 8)), path)


def _cgroup_limits(proc):
    """Read v1/v2 constraints at every ancestor up to the mounted hierarchy root."""
    memberships = []
    for line in _read(proc / "self/cgroup").splitlines():
        _, controllers, path = line.split(":", 2)
        memberships.append((set(controllers.split(",")) - {""}, path))
    observations = []
    visited = set()
    for line in _read(proc / "self/mountinfo").splitlines():
        left, separator, right = line.partition(" - ")
        if not separator:
            continue
        parts, filesystem = left.split(), right.split()
        if filesystem[0] not in {"cgroup", "cgroup2"}:
            continue
        mount_root, mount = _unescape_mount(parts[3]), Path(_unescape_mount(parts[4]))
        controllers = set(filesystem[2].split(","))
        for member_controllers, member_path in memberships:
            v2 = filesystem[0] == "cgroup2"
            if (v2 and member_controllers) or (not v2 and not controllers & member_controllers):
                continue
            if member_path == "/":  # A private cgroup namespace exposes its root as /.
                relative = ""
            elif mount_root == "/":
                relative = member_path.lstrip("/")
            elif member_path == mount_root or member_path.startswith(mount_root.rstrip("/") + "/"):
                relative = member_path[len(mount_root) :].lstrip("/")
            else:
                continue
            current = mount / relative
            for root in (current, *current.parents):
                if root != mount and mount not in root.parents:
                    break
                if root in visited:
                    continue
                visited.add(root)
                limit = _integer(_read(root / ("memory.max" if v2 else "memory.limit_in_bytes")))
                used_text = _read(root / ("memory.current" if v2 else "memory.usage_in_bytes"))
                available = None
                if limit is not None and used_text.isdigit():
                    stats = dict(line.split() for line in _read(root / "memory.stat").splitlines())
                    inactive = int(stats.get("inactive_file" if v2 else "total_inactive_file", 0))
                    available = max(0, limit - max(0, int(used_text) - inactive))
                if v2:
                    cpu = _read(root / "cpu.max").split()
                else:
                    cpu = [_read(root / "cpu.cfs_quota_us"), _read(root / "cpu.cfs_period_us")]
                quota = None
                if len(cpu) == 2 and _integer(cpu[0]) and _integer(cpu[1]):
                    quota = int(cpu[0]) / int(cpu[1])
                if limit is not None or quota is not None:
                    observations.append(
                        {
                            "path": str(root),
                            "memory_limit_bytes": limit,
                            "memory_available_bytes": available,
                            "cpu_quota": quota,
                        }
                    )
    return observations


def detect_capacity(*, proc=Path("/proc"), environ=None):
    """Observe capacity available to this process, not the cluster node's full RAM."""
    environ = os.environ if environ is None else environ
    total, available = _host_memory(proc)
    cpus = os.cpu_count() or 1
    affinity = None
    if hasattr(os, "sched_getaffinity"):
        affinity = len(os.sched_getaffinity(0))
        cpus = min(cpus, affinity)
    slurm_cpus = _integer(environ.get("SLURM_CPUS_PER_TASK"))
    if slurm_cpus:
        cpus = min(cpus, slurm_cpus)
    slurm_memory = _integer(environ.get("SLURM_MEM_PER_NODE"))
    if slurm_memory is None and slurm_cpus:
        per_cpu = _integer(environ.get("SLURM_MEM_PER_CPU"))
        if per_cpu:
            slurm_memory = per_cpu * slurm_cpus
    memory_limits = [total]
    if slurm_memory:
        memory_limits.append(slurm_memory * MIB)
    observations = _cgroup_limits(proc) if sys.platform.startswith("linux") else []
    for observation in observations:
        if observation["memory_limit_bytes"] is not None:
            memory_limits.append(observation["memory_limit_bytes"])
        if observation["memory_available_bytes"] is not None:
            available = min(available, observation["memory_available_bytes"])
        if observation["cpu_quota"] is not None:
            cpus = min(cpus, max(1, math.floor(observation["cpu_quota"])))
    capacity = min(memory_limits)
    return {
        "memory_capacity_bytes": capacity,
        "memory_available_bytes": min(available, capacity),
        "cpu_slots": max(1, cpus),
        "host_memory_bytes": total,
        "affinity_cpus": affinity,
        "slurm_cpus_per_task": slurm_cpus,
        "slurm_memory_mib": slurm_memory,
        "cgroup_observations": observations,
        "availability_note": "Point-in-time estimate; inactive file cache can be reclaimed",
    }


def plan_laptop(
    *,
    memory_fraction=0.35,
    cpu_fraction=0.5,
    memory_gib=None,
    threads=None,
    chunk_mib=None,
    capacity=None,
):
    """Choose shared-process threads and reference pieces; never infer jobs per GB."""
    for name, fraction in (("memory", memory_fraction), ("CPU", cpu_fraction)):
        if not math.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError(f"{name} fraction must be in (0, 1]")
    for name, value in (("threads", threads), ("chunk MiB", chunk_mib)):
        if value is not None and (
            not isinstance(value, int) or isinstance(value, bool) or value < 1
        ):
            raise ValueError(f"{name} must be a positive integer")
    if memory_gib is not None and (not math.isfinite(memory_gib) or memory_gib <= 0):
        raise ValueError("Memory budget must be finite and positive")
    capacity = detect_capacity() if capacity is None else capacity
    available = capacity["memory_available_bytes"]
    ceiling = min(capacity["memory_capacity_bytes"], int(available * 0.8))
    budget = (
        int(memory_gib * GIB)
        if memory_gib is not None
        else min(
            int(capacity["memory_capacity_bytes"] * memory_fraction),
            ceiling,
        )
    )
    if budget > ceiling:
        raise ValueError("Memory budget exceeds available headroom; close other work or lower it")
    if budget < 128 * MIB:
        raise ValueError("Insufficient available memory for laptop planning")
    slots = capacity["cpu_slots"]
    if threads is not None and threads > slots:
        raise ValueError("Requested threads exceed the CPU slots available to this process")
    threads = threads or max(1, math.floor(slots * cpu_fraction))
    chunk_ceiling = max(1, budget // (16 * MIB))
    if chunk_mib is not None and chunk_mib > chunk_ceiling:
        raise ValueError("Reference piece is too large for the planning budget; lower chunk MiB")
    return {
        "schema": "fastalake.laptop-plan.v1",
        "capacity": capacity,
        "memory_budget_bytes": budget,
        "memory_fraction": memory_fraction,
        "cpu_fraction": cpu_fraction,
        "threads_per_process": threads,
        "parallel_processes": 1,
        "reference_chunk_mib": chunk_mib or min(256, chunk_ceiling),
        "blas_threads": 1,
        "scope": (
            "Planning target, not an enforced RAM cap. Steps and acquisitions stay sequential; "
            "Rust threads share process data."
        ),
        "enforcement": (
            "For a hard job limit use Docker --memory or Slurm --mem. "
            "Chunk size does not bound razor, Sage or quantification RAM."
        ),
    }
