"""Record command resource use without adding a monitoring dependency."""

import os
import subprocess
import sys
import time
from pathlib import Path


def local_interpreter(path):
    """Keep a Python environment's launch path, including its executable symlink.

    Resolving ``venv/bin/python`` to the base interpreter discards the venv's
    package environment. See https://docs.python.org/3/library/venv.html.
    Validate the executable while preserving the path the user selected.
    """
    executable = Path(path).expanduser().absolute()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"Python interpreter is not an executable file: {path}")
    return executable


def run_recorded(command, **kwargs):
    """Record POSIX child usage when available.

    ru_maxrss is the largest process RSS reported for the reaped command and its
    waited-for descendants, not simultaneous memory summed across a process tree.
    A Slurm/cgroup memory limit is the separate end-to-end allocation test.
    """
    started = time.monotonic()
    usage = None
    with subprocess.Popen(command, **kwargs) as process:
        try:
            if hasattr(os, "wait4"):
                _, status, usage = os.wait4(process.pid, 0)
                process.returncode = os.waitstatus_to_exitcode(status)
            else:
                process.wait()
        except BaseException:
            process.kill()
            process.wait()
            raise
    rss = None
    if usage is not None and sys.platform in {"linux", "darwin"}:
        rss = int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024))
    return {
        "exit_code": process.returncode,
        "seconds": time.monotonic() - started,
        "max_process_rss_bytes": rss,
        "user_cpu_seconds": usage.ru_utime if usage else None,
        "system_cpu_seconds": usage.ru_stime if usage else None,
        "rss_scope": "Largest process RSS; not summed simultaneous process-tree memory",
    }
