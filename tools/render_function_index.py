#!/usr/bin/env python3
"""Generate source-linked Python and Rust function indexes without importing the package."""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "https://github.com/MannLabs/fasta-lake/blob/main/"
BT = chr(96)
GROUPS = [
    (
        "1. Prepare inputs and plan resources",
        {
            "downloads",
            "capacity",
            "lake_builder",
            "lake_merge",
            "hashing",
            "helpers",
            "dianovo",
            "headers",
            "canonical",
            "contaminants",
            "complexity",
        },
    ),
    ("2. Select acquisition databases", {"chunked", "inference", "ics", "blosum", "refine"}),
    ("3. Add molecular evidence", {"multi_omics"}),
    ("4. Search, group and quantify", {"search", "study", "quantification", "aggregate"}),
    (
        "5. Inspect and interpret results",
        {
            "downstream",
            "qc_plots",
            "quality",
            "taxonomy",
            "taxonomy_reference",
            "annotation",
            "embeddings",
            "networks",
            "network_annotations",
            "community_plots",
            "metaproteomics",
            "exploration",
            "proteomics",
        },
    ),
    (
        "6. Run the workflow and configure tools",
        {
            "tools",
            "cli",
            "run_config",
            "config",
            "presets",
            "resources",
            "rust",
            "provenance",
            "gzip_io",
            "annotation_cli",
            "download_cli",
            "exploration_cli",
            "network_cli",
            "taxonomy_cli",
        },
    ),
    ("7. Diagnostics and comparisons", set()),
]


def functions(tree, prefix=""):
    """Yield qualified names and AST nodes, including nested helpers and methods."""
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            name = prefix + node.name
            if not isinstance(node, ast.ClassDef):
                yield name, node
            yield from functions(node, name + ".")
        else:
            yield from functions(node, prefix)


def cell(text):
    """Keep source summaries on one Markdown table line."""
    return text.replace("|", r"\|").replace("\n", " ").strip()


def python_index():
    """Render package functions and the two main workflow tools in a stable order."""
    files = sorted((ROOT / "fasta_lake").rglob("*.py"))
    files += [ROOT / "tools" / name for name in ("run_manifest.py", "group_study.py")]
    buckets = [[] for _ in GROUPS]
    count, undocumented = 0, []
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        nodes = list(functions(ast.parse(path.read_text())))
        if not nodes:
            continue
        module = relative[:-3].replace("/", ".")
        key = "tools" if relative.startswith("tools/") else module.split(".")[1]
        group = next((i for i, (_, names) in enumerate(GROUPS) if key in names), len(GROUPS) - 1)
        entries = []
        for name, node in sorted(nodes, key=lambda pair: (pair[0].casefold(), pair[1].lineno)):
            doc = ast.get_docstring(node)
            if not doc:
                undocumented.append(f"{relative}:{node.lineno} {name}")
            summary = doc.splitlines()[0] if doc else "See source; documentation pending."
            entries.append(
                f"| [{BT}{name}{BT}]({SOURCE}{relative}#L{node.lineno}) | {cell(summary)} |"
            )
        count += len(entries)
        buckets[group].append((module, entries))
    if undocumented:
        raise ValueError("Missing function docstrings:\n" + "\n".join(undocumented))
    lines = [
        "# Python function index",
        "",
        "For a first reading, use [functions in execution order](function-map.md).",
        f"This generated index contains **{count} functions and methods**, including private",
        "and nested helpers, from the Python package and the two main workflow tools.",
        "Sections follow workflow purpose; modules and names are alphabetical within them.",
        "Links open the exact definition. The [Python API](api/index.md) renders full",
        "docstrings. Legacy and experimental functions remain listed; their presence",
        "here does not make them part of the default runner.",
        "",
        f"Regenerate with {BT}python tools/render_function_index.py --write{BT}.",
        "CI checks the index against the source.",
        "",
    ]
    for (title, _), modules in zip(GROUPS, buckets):
        lines += [f"## {title}", ""]
        for module, entries in sorted(modules):
            lines += [f"### {module}", "", "| Function | Purpose |", "|---|---|", *entries, ""]
    return "\n".join(lines).rstrip() + "\n"


def rust_index():
    """Index tracked Rust declarations outside conventional in-file test modules."""
    tracked = subprocess.check_output(["git", "ls-files", "rust"], cwd=ROOT, text=True).splitlines()
    groups, count = [], 0
    pattern = re.compile(r"^\s*(?:pub(?:\([^)]*\))?\s+)?(?:async\s+)?fn\s+(\w+)")
    for relative in sorted(tracked):
        if not relative.endswith(".rs") or "/tests/" in relative:
            continue
        source = (ROOT / relative).read_text()
        source = re.split(r"#\[cfg\(test\)\]\s*(?:pub\s+)?mod\s+", source, maxsplit=1)[0]
        rows = []
        for number, line in enumerate(source.splitlines(), 1):
            match = pattern.match(line)
            if match:
                rows.append((match.group(1), number))
        if rows:
            groups.append((relative, sorted(rows)))
            count += len(rows)
    lines = [
        "# Rust function index",
        "",
        "Use [functions in execution order](function-map.md) to find the engine for a step.",
        f"This index links **{count} function declarations**, sorted by file and name.",
        "It omits test directories and conventional in-file test modules.",
        "It is a source navigation aid, not a list of public or enabled APIs.",
        "Repeated names belong to different files or implementations.",
        "",
        f"Regenerate with {BT}python tools/render_function_index.py --write{BT}.",
        "",
    ]
    for file, rows in groups:
        lines += [f"## {file}", "", "| Function | Definition |", "|---|---|"]
        lines += [
            f"| {BT}{name}{BT} | [Line {number}]({SOURCE}{file}#L{number}) |"
            for name, number in rows
        ]
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main():
    """Write both indexes, or fail when checked-in copies differ from the source."""
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = {"function-index.md": python_index(), "rust-function-index.md": rust_index()}
    stale = []
    for name, text in expected.items():
        path = ROOT / "docs" / "reference" / name
        if args.write:
            path.write_text(text)
        elif not path.exists() or path.read_text() != text:
            stale.append(name)
    if stale:
        print("Stale function indexes: " + ", ".join(stale))
        print("Run: python tools/render_function_index.py --write")
        return 1
    print("Function indexes written." if args.write else "Function indexes match source.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
