"""Every module shipped in the wheel must import with the core [dev] dependencies.

A module that raises on import is dead on arrival for anyone who reaches it, and
no other test notices unless it imports that exact module. headers/report.py
shipped for weeks importing a module that did not exist. This test walks the
whole package so that cannot happen again.
"""

from __future__ import annotations

import importlib
import pkgutil

import pytest

import fasta_lake

MODULES = sorted(m.name for m in pkgutil.walk_packages(fasta_lake.__path__, "fasta_lake."))


def test_package_has_modules():
    assert len(MODULES) > 50, MODULES


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    importlib.import_module(name)
