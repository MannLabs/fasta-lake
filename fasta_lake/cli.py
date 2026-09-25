"""
Click CLI for FASTA Lake.

Provides command groups for each pipeline stage:

.. code-block:: bash

    fasta-lake lake-build -c fastalake.json
    fasta-lake evidence-lake -d lake.fasta -p predictions/ -o evidence.fasta
    fasta-lake extract-sample -d evidence.fasta -p sample.csv -o sample.fasta
    fasta-lake infer --sample SAMPLE_001 --stage2-fasta s2.fasta \\
        --peptides sample.csv -o output/ --strategy species_budget
    fasta-lake validate fastalake.json
    fasta-lake init --template full
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from fasta_lake import __version__
from fasta_lake.config import create_example_config, load_config
from fasta_lake.entrapment_stages import VALID_ENTRAPMENT_STAGES
from fasta_lake.inference.runner import run_inference
from fasta_lake.inference.strategies import STRATEGIES
from fasta_lake.presets import PRESETS

logger = logging.getLogger(__name__)


def _setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI usage."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )


@click.group()
@click.version_option(version=__version__, prog_name="fasta-lake")
def cli() -> None:
    """FASTA Lake: Sample-specific protein database construction for metaproteomics."""


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------


@cli.group("preset")
def preset_grp() -> None:
    """Inspect or export validated parameter presets (v5).

    Presets are named bundles of stage-specific parameters validated against
    entrapment FDP on 592 MicrobPredict samples. See DESIGN_V5.md.
    """


@preset_grp.command("list")
def preset_list() -> None:
    """List registered parameter presets and their search thresholds."""
    click.echo(f"{'PRESET':<12s} {'search q':>10s}  DESCRIPTION")
    click.echo("-" * 90)
    for name, p in PRESETS.items():
        click.echo(f"{name:<12s} {p.search_fdr:>10g}  {p.description}")
    click.echo("\nThresholds are preset parameters, not whole-workflow error calibration.")


@preset_grp.command("show")
@click.argument("name", type=click.Choice(list(PRESETS)))
def preset_show(name: str) -> None:
    """Show all parameters for a preset."""
    from fasta_lake.presets import describe_preset

    click.echo(describe_preset(name))


@preset_grp.command("sage-config")
@click.argument("name", type=click.Choice(list(PRESETS)))
@click.option(
    "-o", "--output", type=click.Path(), default=None, help="Output JSON path (default: stdout)"
)
@click.option(
    "--base",
    type=click.Path(exists=True),
    default=None,
    help="Base SAGE config to modify (default: bundled sage_dda.json)",
)
def preset_sage(name: str, output: str | None, base: str | None) -> None:
    """Emit a SAGE config with a preset applied."""
    import json

    from fasta_lake.presets import preset_sage_config

    if base is None:
        base = str(Path(__file__).parent / "configs" / "sage_dda.json")
    with open(base) as f:
        base_cfg = json.load(f)

    cfg = preset_sage_config(base_cfg, name)
    out_json = json.dumps(cfg, indent=2)
    if output:
        Path(output).write_text(out_json)
        click.echo(f"Wrote {output}")
    else:
        click.echo(out_json)


# ---------------------------------------------------------------------------
# Contaminants
# ---------------------------------------------------------------------------


@cli.group("contaminants")
def contaminants_grp() -> None:
    """Inspect or apply the contaminant database (v5.1).

    FastaLake bundles a minimal contaminant accession list (keratins,
    trypsin, BSA, caseins, IgGs, etc.). Proteins matching this list are
    tagged with a `CON_` prefix; downstream FDR accounting can exclude
    them via `fasta_lake.contaminants.exclude_contaminants()`.
    """


@contaminants_grp.command("list")
def contam_list() -> None:
    """List the bundled contaminant accessions."""
    from fasta_lake.contaminants import DEFAULT_CONTAMINANTS

    click.echo(f"{len(DEFAULT_CONTAMINANTS)} bundled contaminants:")
    for acc in sorted(DEFAULT_CONTAMINANTS):
        click.echo(f"  {acc}")


@contaminants_grp.command("tag")
@click.option(
    "-i",
    "--input",
    "input_fasta",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTA to scan and tag.",
)
@click.option(
    "-o",
    "--output",
    "output_fasta",
    required=True,
    type=click.Path(),
    help="New output FASTA with CON_-prefixed headers; existing files are refused.",
)
@click.option(
    "--contaminant-db",
    default=None,
    type=click.Path(exists=True),
    help="Optional extra contaminant list (text: 1 accession/line, or FASTA).",
)
@click.option("-v", "--verbose", is_flag=True)
def contam_tag(
    input_fasta: str, output_fasta: str, contaminant_db: str | None, verbose: bool
) -> None:
    """Scan a FASTA, prefix CON_ on headers matching the contaminant list."""
    _setup_logging(verbose)
    from fasta_lake.contaminants import ContaminantDB

    db = ContaminantDB()
    if contaminant_db:
        path = Path(contaminant_db)
        if path.suffix.lower() in (".fasta", ".fa", ".faa"):
            db.load_from_fasta(path)
        else:
            db.load(path)
        click.echo(
            f"Extended contaminant DB with entries from {path} (now {len(db.accessions)} total)"
        )

    n_in = 0
    n_tagged = 0
    with open(input_fasta) as fin, open(output_fasta, "x") as fout:
        for line in fin:
            if line.startswith(">"):
                n_in += 1
                header = line[1:].rstrip("\n")
                pid = header.split()[0]
                if db.is_contaminant(pid) and not pid.startswith(db.prefix):
                    new_pid = f"{db.prefix}{pid}"
                    header = new_pid + header[len(pid) :]
                    n_tagged += 1
                fout.write(">" + header + "\n")
            else:
                fout.write(line)

    click.echo(
        f"Scanned {n_in} proteins. Tagged {n_tagged} as contaminants. Output: {output_fasta}"
    )


# ---------------------------------------------------------------------------
# Entrapment — generate, inject, measure FDP
# ---------------------------------------------------------------------------


@cli.group("entrapment")
def entrapment_grp() -> None:
    """Generate, inject and measure entrapment sequences.

    Entrapment gives an FDR estimate independent of target-decoy: sequences
    known to be absent from the sample are put into the database, and anything
    that matches them is a false discovery.

    \b
    WHERE entrapment is injected decides WHAT is being measured. Four points
    exist, earliest (hardest to survive) first:
      prelake    into the input lake BEFORE the cohort evidence filter. Faces
                 every curation step: the fullest construction test.
      prebuild   onto the clustered lake, then RE-CLUSTERED at 97%, so
                 entrapment must also survive clustering.
      preselect  onto the clustered lake, no re-clustering (DEFAULT). Must be
                 nominated by the sample's own de novo peptides.
      appended   onto the finished per-sample database AFTER parsimony,
                 entering untested. Measures SEARCH error, not construction.

    prebuild and preselect are NOT the same stage under two names. All are
    legitimate, none are the same quantity, so these commands record the stage in
    every output and refuse to pool across stages or classes.

    \b
    PREFIXES are not uniform in this repo: ENTRAP_PFUR_/ENTRAP_SHUF_/..., but
    also ENTRAPHUMAN_ (SIHUMIx) and ENTRAPNN_ (near-neighbour), which do NOT
    match ENTRAP_. A wrong pattern scores entrapment as target and drives the FDP
    to zero silently. The default pattern is the common denominator ENTRAP;
    check any lake with `fasta-lake entrapment prefixes <fasta>`.

    \b
    Typical run:
      fasta-lake entrapment classes
      fasta-lake entrapment generate --class shuffled --source lake.fasta \\
          -o entrap.fasta
      fasta-lake entrapment inject -d lake.fasta -e entrap.fasta \\
          -o augmented.fasta --class shuffled
      # ... run the pipeline on augmented.fasta, then search ...
      fasta-lake entrapment fdp --results-dir results/ --target-fasta lake.fasta \\
          -o fdp/
    """


@entrapment_grp.command("classes")
@click.option("--json", "as_json", is_flag=True, help="Emit the full registry as JSON.")
def entrapment_classes_cmd(as_json: bool) -> None:
    """List entrapment classes with their difficulty and interpretation.

    The classes do NOT measure the same thing. Read the caveats before comparing
    any two numbers produced with different classes.
    """
    import json as _json

    from fasta_lake.entrapment_classes import ENTRAPMENT_CLASSES, class_names

    if as_json:
        click.echo(
            _json.dumps({n: ENTRAPMENT_CLASSES[n].as_metadata() for n in class_names()}, indent=2)
        )
        return

    for name in class_names():
        c = ENTRAPMENT_CLASSES[name]
        flag = "  [recommended headline]" if c.recommended_headline else ""
        click.echo(f"\n{name}{flag}")
        click.echo(f"  prefix:     {c.prefix}")
        click.echo(f"  difficulty: {c.difficulty}")
        click.echo(f"  estimates:  {c.estimates}")
        click.echo(f"  what:       {c.description}")
        click.echo(f"  CAVEAT:     {c.caveat}")
        if c.default_source:
            click.echo(f"  source:     {c.default_source}")
        if c.requires_source:
            click.echo("  source:     REQUIRED (--source)")
        if not c.poolable:
            click.echo("  pooling:    NOT poolable with the foreign-organism classes")
    click.echo()


@entrapment_grp.command("prefixes")
@click.argument("fasta", type=click.Path(exists=True))
@click.option(
    "--candidate", default="ENTRAP", help="Substring identifying an entrapment tag. Default ENTRAP."
)
@click.option("--json", "as_json", is_flag=True, help="Emit as JSON.")
def entrapment_prefixes(fasta: str, candidate: str, as_json: bool) -> None:
    """Report the entrapment prefixes actually present in a FASTA.

    Use this before computing an FDP on a database you did not build. A matcher
    configured for the wrong convention scores every entrapment record as target
    and drives the FDP toward zero, silently, with the run completing normally.
    """
    import json as _json

    from fasta_lake.entrapment_workflow import detect_prefixes

    d = detect_prefixes(fasta, candidate=candidate)
    if as_json:
        click.echo(_json.dumps(d, indent=2))
        return

    click.echo(f"File:     {d['path']}")
    click.echo(f"Headers:  {d['n_headers']:,}")
    click.echo(f"Matching {d['candidate']!r}: {d['n_matching_candidate']:,}")
    if not d["observed"]:
        click.echo(
            "\nNo entrapment-like prefixes found. If this database is "
            "supposed to carry entrapment, the convention differs from "
            f"{candidate!r} -- rerun with --candidate."
        )
        return
    click.echo("\nObserved prefixes:")
    for p, n in d["observed"].items():
        tag = "known" if p in d["known"] else "UNKNOWN convention"
        click.echo(f"  {p:<22s} {n:>10,}  [{tag}]")
    if d["unknown"]:
        click.echo(
            "\nWARNING: unrecognised conventions present: "
            f"{', '.join(d['unknown'])}. Pass --entrap-prefix explicitly "
            "to `entrapment fdp` so they are counted as entrapment."
        )


@entrapment_grp.command("generate")
@click.option(
    "--class",
    "entrap_class",
    default="shuffled",
    help="Entrapment class (default: shuffled, the recommended headline "
    "class). See `fasta-lake entrapment classes`.",
)
@click.option(
    "-i",
    "--source",
    "source_fasta",
    default=None,
    type=click.Path(exists=True),
    help="Source FASTA. Required for --class shuffled (the sequences to "
    "shuffle -- normally the lake you will inject into) and for "
    "classes with no bundled default. Overrides the class default.",
)
@click.option(
    "-o",
    "--output",
    "output_fasta",
    required=True,
    type=click.Path(),
    help="Output entrapment FASTA.",
)
@click.option(
    "--entrapment-seed",
    "seed",
    default=42,
    type=int,
    help="Random seed; same seed gives a byte-identical file. Default 42. Zero is rejected.",
)
@click.option(
    "-n",
    "--n-proteins",
    default=None,
    type=int,
    help="Write only N proteins, chosen by seeded subsample. Default: all.",
)
@click.option(
    "--prefix", default=None, help="Override the class header prefix (default: the class's own)."
)
@click.option(
    "--preserve-cleavage-sites/--no-preserve-cleavage-sites",
    default=True,
    help="Hold K/R in place when shuffling so tryptic peptide lengths "
    "match the source. Default: on.",
)
@click.option("-v", "--verbose", is_flag=True)
def entrapment_generate(
    entrap_class: str,
    source_fasta: str | None,
    output_fasta: str,
    seed: int,
    n_proteins: int | None,
    prefix: str | None,
    preserve_cleavage_sites: bool,
    verbose: bool,
) -> None:
    """Produce an entrapment FASTA, deterministically.

    Headers are tagged with the class prefix (ENTRAP_SHUF_, ENTRAP_PFUR_, ...) so
    downstream hits can be attributed. The same seed and input always produce a
    byte-identical file.
    """
    _setup_logging(verbose)
    from fasta_lake.entrapment_classes import get_class
    from fasta_lake.entrapment_workflow import generate_entrapment_fasta

    cls = get_class(entrap_class)
    stats = generate_entrapment_fasta(
        output_fasta,
        entrap_class=entrap_class,
        source_fasta=source_fasta,
        seed=seed,
        n_proteins=n_proteins,
        prefix=prefix,
        preserve_cleavage_sites=preserve_cleavage_sites,
    )
    click.echo(f"Class:      {stats.entrap_class} ({cls.difficulty})")
    click.echo(f"Estimates:  {cls.estimates}")
    click.echo(f"Source:     {stats.source} ({stats.n_source:,} proteins)")
    click.echo(f"Written:    {stats.n_written:,} sequences, prefix {stats.prefix}")
    click.echo(f"Seed:       {stats.seed}")
    click.echo(f"Output:     {stats.output}")
    click.echo(f"\nCAVEAT: {cls.caveat}")


@entrapment_grp.command("inject")
@click.option(
    "-d",
    "--database",
    required=True,
    type=click.Path(exists=True),
    help="Database to inject into. For the pre-selection stages this is "
    "the lake that de novo selection will run against; for appended "
    "it is the finished per-sample parsimony FASTA.",
)
@click.option(
    "-e",
    "--entrapment",
    required=True,
    type=click.Path(exists=True),
    help="Entrapment FASTA from `entrapment generate`.",
)
@click.option(
    "-o",
    "--output",
    "output_fasta",
    required=True,
    type=click.Path(),
    help="Augmented FASTA to write.",
)
@click.option(
    "--stage",
    type=click.Choice(list(VALID_ENTRAPMENT_STAGES)),
    default="preselect",
    help="Injection stage, earliest first: prelake (before the evidence "
    "filter), prebuild (re-clustered with the lake), preselect "
    "(default; clustered lake, no re-clustering), appended (after "
    "parsimony). The first three measure database-construction "
    "error and differ by how much curation entrapment survived; "
    "appended measures search error. Recorded in the manifest so "
    "the result cannot be misread later.",
)
@click.option(
    "--class",
    "entrap_class",
    default=None,
    help="Entrapment class, recorded in the manifest along with its difficulty and interpretation.",
)
@click.option(
    "--entrapment-seed",
    "seed",
    default=None,
    type=int,
    help="Seed used to generate the entrapment, for the manifest.",
)
@click.option(
    "--entrap-prefix",
    "entrap_prefix",
    default=None,
    help="Header prefix identifying entrapment. Defaults to the class "
    "prefix, or ENTRAP when no class is given. Recorded in the "
    "manifest and reused by `entrapment fdp`.",
)
@click.option(
    "--manifest",
    "manifest_path",
    default=None,
    type=click.Path(),
    help="Manifest path (default: entrapment_manifest.json beside output).",
)
@click.option("--no-checksums", is_flag=True, help="Skip SHA-256 of the inputs.")
@click.option("-v", "--verbose", is_flag=True)
def entrapment_inject(
    database: str,
    entrapment: str,
    output_fasta: str,
    stage: str,
    entrap_class: str | None,
    seed: int | None,
    entrap_prefix: str | None,
    manifest_path: str | None,
    no_checksums: bool,
    verbose: bool,
) -> None:
    """Concatenate entrapment into a database and write a manifest.

    The manifest records the stage, the class, the seed, the entrapment:target
    ratio and the counts. Protein counts come from header lines only.
    """
    _setup_logging(verbose)
    from fasta_lake.entrapment_workflow import inject_entrapment

    m = inject_entrapment(
        database,
        entrapment,
        output_fasta,
        stage=stage,
        entrap_class=entrap_class,
        seed=seed,
        prefix=entrap_prefix,
        manifest_path=manifest_path,
        checksums=not no_checksums,
    )
    click.echo(f"Stage:      {m['stage']}  ({m['measures']})")
    click.echo(f"Class:      {m['class']}")
    click.echo(f"Target:     {m['n_target']:,} proteins")
    click.echo(f"Entrapment: {m['n_entrapment']:,} proteins")
    click.echo(f"Ratio E:T:  {m['ratio_entrapment_to_target']}")
    click.echo(f"Prefix:     {m['entrapment_match_pattern']}")
    click.echo(f"Output:     {m['output_fasta']}")
    click.echo(f"Manifest:   {m['manifest_path']}")
    click.echo(f"\n{m['stage_meaning']}")


@entrapment_grp.command("fdp")
@click.option(
    "--results-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory of per-sample search results "
    "(SAMPLE/results.sage.tsv or SAMPLE/*/results.sage.tsv).",
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Output directory for the per-sample TSV and the FDP record.",
)
@click.option("--dataset", default="", help="Dataset/cohort name for the record.")
@click.option(
    "--class",
    "entrap_class",
    default=None,
    help="Entrapment class. Read from the manifest when present; an "
    "explicit value that contradicts the manifest is refused.",
)
@click.option(
    "--stage",
    type=click.Choice(list(VALID_ENTRAPMENT_STAGES)),
    default=None,
    help="Injection stage. Read from the manifest when present. Required "
    "if no manifest exists -- it will not be guessed.",
)
@click.option(
    "--target-fasta",
    default=None,
    type=click.Path(exists=True),
    help="Target sequence space for the conserved-peptide test. Every "
    "entrapment-attributed peptide is checked for occurrence here "
    "(I/L normalised, exact substring); those found are genuinely "
    "conserved, not false discoveries. Required for homologous "
    "classes such as plant.",
)
@click.option(
    "--no-conserved-check",
    is_flag=True,
    help="Skip the conserved-peptide correction and report an explicitly UNCORRECTED FDP.",
)
@click.option(
    "--normalize-il/--no-normalize-il",
    default=True,
    help="Treat I and L as equivalent when matching. Default: on.",
)
@click.option("--q-threshold", default=0.01, type=float, help="FDR cutoff. Default 0.01.")
@click.option("--q-column", default="peptide_q", help="FDR column. Default peptide_q.")
@click.option("--min-length", default=9, type=int, help="Peptide length floor. Default 9.")
@click.option(
    "--entrap-prefix",
    "entrap_prefix",
    default=None,
    help="Header prefix identifying entrapment. Defaults to what the "
    "manifest recorded, else the class prefix. Set this when the "
    "run used a non-standard convention such as ENTRAPHUMAN_ or "
    "ENTRAPNN_ -- a wrong pattern silently reports ~0% FDP.",
)
@click.option(
    "--decoy-prefix",
    default="rev_",
    help="Decoy tag, used only if the results have no `label` column. SAGE writes rev_.",
)
@click.option(
    "--db-pattern",
    default=None,
    help="Glob (relative to a sample dir) for that sample's searched "
    "FASTA, e.g. '*_razor.fasta'. Enables per-sample database "
    "partition sizes and the size-corrected estimators.",
)
@click.option(
    "--results-name",
    default="results.sage.tsv",
    help="Results filename, or '<subdir>/results.sage.tsv' to choose one of several "
    "entrapment runs inside a sample directory.",
)
@click.option("--peptide-level", is_flag=True, help="Count distinct peptides instead of spectra.")
@click.option(
    "--allow-mixed",
    is_flag=True,
    help="Permit pooling across injection stages or classes. Produces a "
    "number that mixes different quantities; use only deliberately.",
)
@click.option(
    "--entrapment-seed", "seed", default=42, type=int, help="Seed for the bootstrap CI. Default 42."
)
@click.option("--bootstrap-draws", default=2000, type=int, help="Bootstrap draws. Default 2000.")
@click.option("-v", "--verbose", is_flag=True)
def entrapment_fdp(
    results_dir: str,
    output_dir: str,
    dataset: str,
    entrap_class: str | None,
    stage: str | None,
    target_fasta: str | None,
    no_conserved_check: bool,
    normalize_il: bool,
    q_threshold: float,
    q_column: str,
    min_length: int,
    entrap_prefix: str | None,
    decoy_prefix: str,
    db_pattern: str | None,
    results_name: str,
    peptide_level: bool,
    allow_mixed: bool,
    seed: int,
    bootstrap_draws: int,
    verbose: bool,
) -> None:
    """Compute the FDP from search results, with n, median, IQR and bootstrap CI.

    Reports raw and conserved-corrected FDP together, and states how many samples
    were evaluable. A sample with no entrapment sequences in its searched database
    has an UNDEFINED FDP, not a zero, and is excluded rather than counted as 0.0.
    """
    _setup_logging(verbose)
    from fasta_lake.entrapment_workflow import compute_fdp, write_fdp_reports

    report = compute_fdp(
        results_dir,
        dataset=dataset,
        entrap_class=entrap_class,
        stage=stage,
        target_fasta=target_fasta,
        q_threshold=q_threshold,
        min_length=min_length,
        q_column=q_column,
        prefix=entrap_prefix,
        decoy_prefix=decoy_prefix,
        db_pattern=db_pattern,
        results_name=results_name,
        seed=seed,
        n_draws=bootstrap_draws,
        peptide_level=peptide_level,
        allow_mixed=allow_mixed,
        conserved_check=not no_conserved_check,
        normalize_il=normalize_il,
    )
    written = write_fdp_reports(report, output_dir)

    cov = report["coverage"]
    meta = report["class_metadata"]
    click.echo(f"Dataset:    {report['dataset'] or '(unnamed)'}")
    click.echo(f"Class:      {report['class']}  ({meta['difficulty']})")
    click.echo(f"Stage:      {report['stage']}  -> measures {report['measures']}")
    click.echo(f"Estimates:  {meta['estimates']}")
    click.echo(f"Counting:   {report['counting']}")
    click.echo()
    click.echo("=== Coverage ===")
    click.echo(f"  samples attempted:              {cov['n_attempted']}")
    click.echo(f"  samples evaluable:              {cov['n_evaluable']}")
    click.echo(f"  empty entrapment partition:     {cov['n_empty_entrapment_partition']}")
    click.echo(f"  unknown entrapment partition:   {cov['n_unknown_entrapment_partition']}")
    if cov["n_empty_entrapment_partition"]:
        click.echo(
            "  NOTE: samples with an empty entrapment partition have an "
            "UNDEFINED FDP, not 0.0. They are excluded from the summary."
        )

    cons = report["conserved_peptides"]
    click.echo()
    click.echo("=== Conserved-peptide correction ===")
    click.echo(f"  status: {cons['status']}")
    if cons["performed"]:
        click.echo(f"  entrapment peptides tested:     {cons['n_entrapment_peptides_tested']}")
        click.echo(
            f"  also present in target space:   {cons['n_conserved']}"
            f" ({(cons['conserved_rate'] or 0):.1%})"
        )
        nc = cons["negative_control"]
        click.echo(
            f"  reversed-peptide control:       {nc['n_found']}/{nc['n_tested']}"
            f" ({(nc['rate'] or 0):.2%}) — should be ~0"
        )
        if nc["warning"]:
            click.echo(f"  WARNING: {nc['warning']}")
        if "false_exculpation_calibration" in cons:
            click.echo(f"  {cons['false_exculpation_calibration']['note']}")

    click.echo()
    click.echo("=== FDP (per-sample median [IQR], bootstrap 95% CI) ===")
    click.echo(f"{'estimator':<36s} {'n':>4s} {'median':>10s} {'IQR':>22s} {'95% CI':>22s}")
    for key, block in report["summary"].items():
        for kind in ("raw", "conserved_corrected"):
            s = block[kind]
            if not s["n"]:
                continue
            label = key if kind == "raw" else f"{key} (cons.corr)"
            ci = (
                "  n/a"
                if s["ci_lo"] is None
                else f"[{s['ci_lo'] * 100:.4f}, {s['ci_hi'] * 100:.4f}]%"
            )
            iqr = f"[{s['q1'] * 100:.4f}, {s['q3'] * 100:.4f}]%"
            click.echo(f"{label:<36s} {s['n']:>4d} {s['median'] * 100:>9.4f}% {iqr:>22s} {ci:>22s}")

    click.echo()
    click.echo(f"CAVEAT ({report['class']}): {meta['caveat']}")
    click.echo()
    for what, path in written.items():
        click.echo(f"{what:<12s} {path}")


@entrapment_grp.command("prepare")
@click.option(
    "-c",
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True),
    help="FastaLake JSON config with an `entrapment` section.",
)
@click.option(
    "-d",
    "--database",
    required=True,
    type=click.Path(exists=True),
    help="Database to inject into (the lake de novo selection runs against).",
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Directory for entrapment.fasta, augmented_database.fasta and the manifest.",
)
@click.option("-v", "--verbose", is_flag=True)
def entrapment_prepare(config_path: str, database: str, output_dir: str, verbose: bool) -> None:
    """Generate + inject entrapment as configured by a fastalake.json.

    This is the config-driven form of `generate` followed by `inject`. It is a
    no-op when `entrapment.enabled` is false.
    """
    _setup_logging(verbose)
    from fasta_lake.entrapment_workflow import prepare_from_config

    cfg = load_config(config_path)
    m = prepare_from_config(cfg, database, output_dir)
    if m is None:
        click.echo("entrapment.enabled is false — nothing to do.")
        click.echo(f"Use the database unchanged: {database}")
        return
    click.echo(f"Class:      {m['class']}")
    click.echo(f"Stage:      {m['stage']}  ({m['measures']})")
    click.echo(f"Target:     {m['n_target']:,} proteins")
    click.echo(f"Entrapment: {m['n_entrapment']:,} proteins")
    click.echo(f"Augmented:  {m['output_fasta']}")
    click.echo(f"Manifest:   {m['manifest_path']}")


# ---------------------------------------------------------------------------
# Tune — post-hoc parameter recommendation from entrapment results
# ---------------------------------------------------------------------------


@cli.command("tune")
@click.option(
    "--results-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory with per-sample entrapment results. Each sample subdir "
    "must contain entrapment/, entrapment_plant/, and/or entrapment_archaea/ "
    "subdirs with results.sage.tsv inside.",
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Output directory for tuning_report.csv and recommended.json.",
)
@click.option(
    "--target-fdp",
    default=2.0,
    type=float,
    help="Max acceptable shuffled FDP (percent). Default 2.0.",
)
@click.option(
    "--target-db-size",
    default=75000,
    type=int,
    help="Target database size used for r = E/T. Default 75000 (typical MP sample).",
)
@click.option("-v", "--verbose", is_flag=True)
def tune_cmd(
    results_dir: str, output_dir: str, target_fdp: float, target_db_size: int, verbose: bool
) -> None:
    """Recommend SAGE parameters by post-hoc sweep over existing entrapment results.

    Takes a directory of SAGE entrapment searches, applies a grid of
    (q_threshold, min_length, min_score) filters, and recommends the filter
    combination that achieves --target-fdp while maximising target PSMs.

    Runs in seconds on 100s of samples — no re-search required.
    """
    from fasta_lake.tune import tune, write_reports

    _setup_logging(verbose)

    report = tune(
        results_dir=Path(results_dir),
        target_fdp_pct=target_fdp,
        target_db_size_fn=lambda _sid: target_db_size,
    )

    out = Path(output_dir)
    write_reports(report, out)

    click.echo(f"Samples analyzed: {report['n_samples']}")
    click.echo(
        f"Filter grid size: {len(report['all_results']) // 3 if report['all_results'] else 0} "
        f"× 3 entrap types = {len(report['all_results'])} evaluations"
    )
    click.echo(f"Verdict: {report['winner_verdict']}")
    if report["winner"]:
        w = report["winner"]
        click.echo()
        click.echo("=== RECOMMENDED ===")
        click.echo(f"  Filter:      {w['filter']}")
        click.echo(f"  q_threshold: {w['q_threshold']}")
        click.echo(f"  min_length:  {w['min_length']}")
        click.echo(f"  min_score:   {w['min_score']}")
        click.echo(f"  Shuffled FDP (median): {w['median_fdp_pct']:.2f}%")
        click.echo(f"  Target PSMs (median):  {w['median_target']:.0f}")
    click.echo()
    click.echo(f"Full report: {out / 'tuning_report.csv'}")
    click.echo(f"Recommended: {out / 'recommended.json'}")


# ---------------------------------------------------------------------------
# Diagnose — FDR calibration curves + outlier detection
# ---------------------------------------------------------------------------


@cli.command("diagnose")
@click.option(
    "--results-dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory with per-sample entrapment results (same layout as `tune`).",
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Output directory for calibration curve + outlier list.",
)
@click.option(
    "--min-length", default=9, type=int, help="Peptide min length (applied post-hoc). Default 9."
)
@click.option(
    "--z-threshold",
    default=2.0,
    type=float,
    help="Z-score threshold for outlier flag. Default 2.0.",
)
@click.option("-v", "--verbose", is_flag=True)
def diagnose_cmd(
    results_dir: str, output_dir: str, min_length: int, z_threshold: float, verbose: bool
) -> None:
    """Produce FDR calibration curve + outlier sample list.

    For each entrapment type (shuffled, plant, archaea), computes empirical
    combined-method FDP at a grid of nominal q-value thresholds. Writes
    calibration_curve.csv. Also flags samples with outlier FDP for review.
    """
    from fasta_lake.diagnose import calibration_curve, outlier_samples, write_report

    _setup_logging(verbose)

    rdir = Path(results_dir)
    cal = calibration_curve(rdir, min_length=min_length)
    outliers = outlier_samples(rdir, fdp_z_threshold=z_threshold, min_length=min_length)
    report = {"calibration": cal, "outliers": outliers}
    write_report(report, Path(output_dir))

    click.echo(f"Samples analyzed: {cal['n_samples']}")
    click.echo(f"Min length filter: {min_length}")
    click.echo()
    click.echo("=== Calibration (shuffled entrapment, PSM- and protein-level FDP) ===")
    click.echo(
        f"{'nominal_q':>10s} {'PSM_FDP':>9s} {'prot_FDP':>10s} {'ratio':>8s} {'n':>4s} "
        f"{'contams':>8s}"
    )
    for p in cal["curve"]:
        if p["entrap_type"] != "shuffled":
            continue
        click.echo(
            f"{p['nominal_q']:>9.4f}  "
            f"{p['median_psm_fdp_pct']:>7.2f}% "
            f"{p['median_prot_fdp_pct']:>8.2f}% "
            f"{p['calibration_ratio_psm']:>6.2f}x "
            f"{p['n_samples']:>4d} "
            f"{p['total_contaminants']:>8d}"
        )
    click.echo()
    click.echo("  ratio = median_PSM_FDP / nominal_q (should approach 1.0 if calibrated)")
    click.echo("  contams = total PSMs matching bundled contaminant accessions (excluded from FDP)")
    click.echo()
    if outliers:
        click.echo(f"=== Outlier samples (z > {z_threshold}) ===")
        for o in outliers[:10]:
            click.echo(
                f"  {o['sample']:<60s} PSM_FDP={o['shuffled_psm_fdp_pct']:.2f}% "
                f"(z={o['z_score']:+.2f}, {o['target_psms']} targets / "
                f"{o['target_proteins']} proteins)"
            )
        if len(outliers) > 10:
            click.echo(f"  ... ({len(outliers) - 10} more)")
    else:
        click.echo("No outliers flagged.")
    click.echo()
    click.echo(f"Full reports: {Path(output_dir)}/{{calibration_curve.csv,outlier_samples.json}}")


# ---------------------------------------------------------------------------
# Stage 3: Protein inference
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--sample", required=True, help="Sample ID")
@click.option(
    "--stage2-fasta",
    required=True,
    type=click.Path(exists=True),
    help="Stage 2 FASTA file",
)
@click.option(
    "--peptides",
    required=True,
    type=click.Path(exists=True),
    help="AlphaNovo predictions CSV",
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Output directory",
)
@click.option(
    "--strategy",
    default="species_budget",
    type=click.Choice(list(STRATEGIES.keys())),
    help=(
        "Experimental Python strategy (default: species_budget). "
        "See tools/run_manifest.py for the validated Rust workflow."
    ),
)
@click.option("--min-peptides", default=2, type=int, help="Min peptides for uniform_2pep")
@click.option(
    "--taxonomy-map",
    default=None,
    type=click.Path(exists=True),
    help="TSV mapping protein_id -> species. species_budget needs a species for every "
    "protein with evidence; ids without one (assembly proteins, smORFs, content-hash "
    "lakes) must be covered here or the strategy refuses to run.",
)
@click.option(
    "--taxonomy-pruning",
    is_flag=True,
    default=False,
    help="Enable frequentist pruning: shared peptides only assigned to abundant species. Requires "
    "--taxonomy-map.",
)
@click.option(
    "--kmer-rescue-db",
    default=None,
    type=click.Path(exists=True),
    help="FASTA database for k-mer rescue of failed peptides. "
    "Peptides that don't exact-match any Stage 2 protein are rescued "
    "via mass k-mer anchoring against this database. Rescued proteins "
    "are injected before inference runs.",
)
@click.option(
    "--kmer-min-votes", default=5, type=int, help="Min k-mer votes for rescue (default: 5)"
)
@click.option(
    "--kmer-threads", default=4, type=int, help="Threads for mass_kmer_anchor (default: 4)"
)
@click.option(
    "--preset",
    default=None,
    type=click.Choice(list(PRESETS)),
    help="Apply a parameter preset (overrides other --min-peptides / --strategy flags). "
    "See DESIGN_V5.md for preset parameters and empirical FDP for each.",
)
@click.option("-v", "--verbose", is_flag=True, help="Verbose logging")
def infer(
    sample: str,
    stage2_fasta: str,
    peptides: str,
    output_dir: str,
    strategy: str,
    min_peptides: int,
    taxonomy_map: str | None,
    taxonomy_pruning: bool,
    kmer_rescue_db: str | None,
    kmer_min_votes: int,
    kmer_threads: int,
    preset: str | None,
    verbose: bool,
) -> None:
    """Experimental Python inference, distinct from the validated Rust workflow."""
    destination = Path(output_dir)
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise click.ClickException("Output must be new or empty; choose a fresh directory")
    click.echo(
        "Experimental Python inference: evidence weighting and tie rules differ from "
        "Rust razor/hash-acc. Use tools/run_manifest.py for the validated workflow.",
        err=True,
    )
    _setup_logging(verbose)

    # Apply preset if given — overrides individual flags
    min_length = 9  # v5 default
    if preset is not None:
        from fasta_lake.presets import get_preset

        p = get_preset(preset)
        click.echo(
            f"[preset={preset}] min_length={p.denovo_min_length} "
            f"strategy={p.infer_strategy} "
            f"min_peptides={p.infer_min_peptides} "
            f"kmer_rescue={p.kmer_rescue}"
        )
        strategy = p.infer_strategy
        min_peptides = p.infer_min_peptides
        min_length = p.denovo_min_length
        if not p.kmer_rescue:
            kmer_rescue_db = None

    stats = run_inference(
        sample_id=sample,
        stage2_fasta=stage2_fasta,
        alphanovo_csv=peptides,
        output_dir=output_dir,
        strategy=strategy,
        min_peptides=min_peptides,
        min_length=min_length,
        taxonomy_map_path=taxonomy_map,
        taxonomy_pruning=taxonomy_pruning,
        kmer_rescue_db=kmer_rescue_db,
        kmer_min_votes=kmer_min_votes,
        kmer_threads=kmer_threads,
    )

    click.echo(f"Selected proteins: {stats['selected_proteins']:,}")
    click.echo(f"Reduction: {stats['reduction_factor']:.2f}x")
    if stats.get("kmer_proteins_new"):
        click.echo(f"K-mer rescued: {stats['kmer_proteins_new']} new proteins injected")
    click.echo(f"Runtime: {stats['runtime_seconds']:.1f}s")


# ---------------------------------------------------------------------------
# Engine-aware FASTA export
# ---------------------------------------------------------------------------


@cli.command("export")
@click.option(
    "-i",
    "--input",
    "input_fasta",
    required=True,
    type=click.Path(exists=True),
    help="Input inference FASTA from FastaLake",
)
@click.option(
    "-o",
    "--output",
    "output_fasta",
    required=True,
    type=click.Path(),
    help="Engine-ready output FASTA",
)
@click.option(
    "-e",
    "--engine",
    required=True,
    type=click.Choice(["sage", "diann", "msfragger"]),
    help="Target search engine",
)
@click.option(
    "--canonical",
    default=None,
    type=click.Path(exists=True),
    help="SwissProt FASTA for gene resolution (DIA-NN)",
)
@click.option("-v", "--verbose", is_flag=True)
def export(
    input_fasta: str, output_fasta: str, engine: str, canonical: str | None, verbose: bool
) -> None:
    """Prepare a FASTA for a specific search engine.

    Each engine has specific requirements:
      sage      — adds GN= tags, no decoys (SAGE generates internal decoys)
      diann     — gene name as first word of description (DIA-NN requirement)
      msfragger — adds GN= tags + reversed decoy sequences (rev_ prefix)
    """
    _setup_logging(verbose)
    from fasta_lake.headers import prepare_fasta_for_engine

    stats = prepare_fasta_for_engine(
        input_fasta,
        output_fasta,
        engine=engine,
        canonical_fasta=canonical,
    )
    click.echo(f"Engine:   {stats.engine}")
    click.echo(f"Targets:  {stats.target_proteins:,}")
    if stats.decoy_proteins:
        click.echo(f"Decoys:   {stats.decoy_proteins:,}")
    click.echo(f"Total:    {stats.total_proteins:,}")
    click.echo(f"Output:   {output_fasta}")


# ---------------------------------------------------------------------------
# Raw data conversion: .raw → .mzML → .hdf → predictions
# ---------------------------------------------------------------------------


@cli.command("convert")
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Input: .raw file, .mzML file, .hdf file, or directory of any",
)
@click.option(
    "-o",
    "--output",
    "output_dir",
    required=True,
    type=click.Path(),
    help="Output directory for converted files",
)
@click.option(
    "--step",
    type=click.Choice(["raw-to-mzml", "mzml-to-hdf", "hdf-to-csv", "all"]),
    default="all",
    help="Which conversion step(s) to run (default: all)",
)
@click.option(
    "--checkpoint",
    type=click.Choice(["old", "rope"]),
    default="rope",
    help="AlphaNovo checkpoint for de novo inference (default: rope)",
)
@click.option(
    "--validate/--no-validate",
    default=True,
    help="Validate HDF files before inference (catches corruption)",
)
@click.option("-v", "--verbose", is_flag=True)
def convert(
    input_path: str, output_dir: str, step: str, checkpoint: str, validate: bool, verbose: bool
) -> None:
    """Convert raw mass spec data to de novo predictions.

    Handles the full chain: .raw → .mzML → .hdf → predictions.csv

    Each step can be run independently or all at once:

    \b
      raw-to-mzml   Thermo .raw → .mzML (ThermoRawFileParser v2)
      mzml-to-hdf   .mzML → .hdf spectral library (AlphaNovo converter)
      hdf-to-csv    .hdf → predictions.csv (AlphaNovo inference)
      all           Run all steps in sequence

    IMPORTANT: AlphaNovo inference uses save_speclib=false to prevent HDF
    corruption. The --validate flag (on by default) checks that input HDFs
    are clean before running inference.

    Examples:

    \b
      # Full pipeline from raw files
      fasta-lake convert -i data/raw/ -o data/predictions/ --step all

    \b
      # Just convert mzML to HDF
      fasta-lake convert -i data/mzml/ -o data/hdf/ --step mzml-to-hdf

    \b
      # Run inference on clean HDFs
      fasta-lake convert -i data/hdf/ -o data/pred/ --step hdf-to-csv --checkpoint rope
    """
    _setup_logging(verbose)
    import subprocess
    from pathlib import Path

    input_path = Path(input_path)
    output_dir = Path(output_dir)
    import os

    configured_scripts = os.environ.get("FASTALAKE_CONVERSION_SCRIPTS")
    if not configured_scripts or not Path(configured_scripts).is_dir():
        raise click.ClickException(
            "Experimental conversion requires FASTALAKE_CONVERSION_SCRIPTS pointing "
            "to your independently configured conversion and AlphaNovo scripts. "
            "The validated workflow starts from predictions CSV and mzML files."
        )
    if output_dir.exists():
        raise click.ClickException("Output exists; use a fresh conversion directory")
    scripts_dir = Path(configured_scripts)
    output_dir.mkdir(parents=True, exist_ok=False)

    steps_to_run = []
    if step == "all":
        # Detect what we have and determine steps
        if input_path.suffix.lower() == ".raw" or (
            input_path.is_dir() and list(input_path.glob("*.raw"))
        ):
            steps_to_run = ["raw-to-mzml", "mzml-to-hdf", "hdf-to-csv"]
        elif input_path.suffix.lower() == ".mzml" or (
            input_path.is_dir() and list(input_path.glob("*.mzML"))
        ):
            steps_to_run = ["mzml-to-hdf", "hdf-to-csv"]
        elif input_path.suffix.lower() == ".hdf" or (
            input_path.is_dir() and list(input_path.glob("*.hdf"))
        ):
            steps_to_run = ["hdf-to-csv"]
        else:
            click.echo("ERROR: Could not detect input type. Use --step to specify.")
            raise SystemExit(1)
    else:
        steps_to_run = [step]

    current_input = str(input_path)

    for s in steps_to_run:
        click.echo(f"\n{'=' * 60}")
        click.echo(f"Step: {s}")
        click.echo(f"{'=' * 60}")

        if s == "raw-to-mzml":
            script = scripts_dir / "convert_raw_to_mzml.sh"
            step_output = str(output_dir / "mzml")
            Path(step_output).mkdir(exist_ok=True)
            result = subprocess.run(
                ["bash", str(script), "-i", current_input, "-o", step_output],
                capture_output=False,
            )
            if result.returncode != 0:
                click.echo(f"ERROR: raw-to-mzml failed (exit {result.returncode})")
                raise SystemExit(1)
            current_input = step_output

        elif s == "mzml-to-hdf":
            script = scripts_dir / "convert_mzml_to_hdf.sh"
            step_output = str(output_dir / "hdf")
            Path(step_output).mkdir(exist_ok=True)
            result = subprocess.run(
                ["bash", str(script), "-i", current_input, "-o", step_output],
                capture_output=False,
            )
            if result.returncode != 0:
                click.echo(f"ERROR: mzml-to-hdf failed (exit {result.returncode})")
                raise SystemExit(1)
            current_input = step_output

        elif s == "hdf-to-csv":
            script = scripts_dir / "run_denovo.sh"
            step_output = str(output_dir / "predictions")
            Path(step_output).mkdir(exist_ok=True)
            args = ["bash", str(script), "-i", current_input, "-o", step_output, "-c", checkpoint]
            if not validate:
                args.append("--no-validate")
            result = subprocess.run(args, capture_output=False)
            if result.returncode != 0:
                click.echo(f"ERROR: hdf-to-csv failed (exit {result.returncode})")
                raise SystemExit(1)

    click.echo(f"\nConversion complete. Output in: {output_dir}")


# ---------------------------------------------------------------------------
# Stages 0-2: Rust binary wrappers
# ---------------------------------------------------------------------------


@cli.command("lake-build")
@click.option("-c", "--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("-v", "--verbose", is_flag=True)
def lake_build(config_path: str, verbose: bool) -> None:
    """Stage 0: Build reference protein lake from databases."""
    _setup_logging(verbose)
    from fasta_lake.rust.binaries import run_lake_builder

    cfg = load_config(config_path)
    sources = [(s.tag, s.priority, s.path) for s in cfg.database.sources]
    output = Path(cfg.output.directory) / cfg.database.output
    output.parent.mkdir(parents=True, exist_ok=True)

    run_lake_builder(
        sources=sources, output=output, threads=cfg.rust.threads, binary_path=cfg.rust.lake_builder
    )
    click.echo(f"Lake built: {output}")


@cli.command("evidence-lake")
@click.option("-d", "--database", required=True, type=click.Path(exists=True))
@click.option("-p", "--predictions", required=True, type=click.Path(exists=True))
@click.option("-o", "--output", required=True, type=click.Path())
@click.option("-s", "--source-tag", default="REF")
@click.option("--normalize-il/--no-normalize-il", default=True)
@click.option("--top-percent", type=float, default=None)
@click.option("--output-hits", type=click.Path(), default=None)
@click.option("-t", "--threads", type=int, default=None)
@click.option("-v", "--verbose", is_flag=True)
def evidence_lake(
    database: str,
    predictions: str,
    output: str,
    source_tag: str,
    normalize_il: bool,
    top_percent: float | None,
    output_hits: str | None,
    threads: int | None,
    verbose: bool,
) -> None:
    """Stage 1: Create evidence lake from peptide predictions."""
    _setup_logging(verbose)
    from fasta_lake.rust.binaries import run_fasta_extractor

    stats = run_fasta_extractor(
        database=database,
        predictions=predictions,
        output=output,
        source_tag=source_tag,
        normalize_il=normalize_il,
        top_percent=top_percent,
        output_hits=output_hits,
        threads=threads,
    )
    if stats:
        click.echo(f"Proteins matched: {stats.get('proteins_matched', '?'):,}")


@cli.command("extract-sample")
@click.option("-d", "--database", required=True, type=click.Path(exists=True))
@click.option("-p", "--predictions", required=True, type=click.Path(exists=True))
@click.option("-o", "--output", required=True, type=click.Path())
@click.option("-s", "--source-tag", default="SAMPLE")
@click.option("--normalize-il/--no-normalize-il", default=True)
@click.option("--output-hits", type=click.Path(), default=None)
@click.option("-t", "--threads", type=int, default=None)
@click.option(
    "-c",
    "--config",
    "config_path",
    default=None,
    type=click.Path(exists=True),
    help="FastaLake JSON config. When `entrapment.enabled` is true and the "
    "stage is pre-selection, entrapment is generated and injected into "
    "--database first, so entrapment competes for de novo nomination "
    "exactly like a real sequence. A manifest is written next to the "
    "augmented lake.",
)
@click.option(
    "--entrapment-dir",
    default=None,
    type=click.Path(),
    help="Where to put the augmented lake and manifest (default: <output>.entrapment/).",
)
@click.option("-v", "--verbose", is_flag=True)
def extract_sample(
    database: str,
    predictions: str,
    output: str,
    source_tag: str,
    normalize_il: bool,
    output_hits: str | None,
    threads: int | None,
    config_path: str | None,
    entrapment_dir: str | None,
    verbose: bool,
) -> None:
    """Stage 2: Extract per-sample protein database.

    With -c/--config and `entrapment.enabled: true`, entrapment is injected into
    the database before extraction (stage `pre-selection`), which is what makes an
    entrapment hit a measurement of database-construction error rather than of
    search error.
    """
    _setup_logging(verbose)
    from fasta_lake.rust.binaries import run_fasta_extractor

    if config_path:
        database = _maybe_inject_entrapment(
            config_path,
            database,
            output,
            entrapment_dir,
        )

    stats = run_fasta_extractor(
        database=database,
        predictions=predictions,
        output=output,
        source_tag=source_tag,
        normalize_il=normalize_il,
        output_hits=output_hits,
        threads=threads,
    )
    if stats:
        click.echo(f"Proteins matched: {stats.get('proteins_matched', '?'):,}")


def _maybe_inject_entrapment(
    config_path: str,
    database: str,
    output: str,
    entrapment_dir: str | None,
) -> str:
    """Apply a config's entrapment section to a database, returning what to search.

    Returns ``database`` unchanged when entrapment is disabled, otherwise the path
    to the augmented database. Refuses to inject at `appended` here: appending
    happens after parsimony, not before extraction.
    """
    from fasta_lake.entrapment_workflow import STAGE_PRESELECT, prepare_from_config

    cfg = load_config(config_path)
    ec = cfg.entrapment
    if not ec.enabled:
        return database

    if ec.stage != STAGE_PRESELECT:
        raise click.ClickException(
            f"entrapment.stage is {ec.stage!r}, which is applied AFTER parsimony, "
            "not before per-sample extraction. Either set the stage to "
            "'preselect' (the default, and the setting that makes an entrapment "
            "hit a database-construction error), or inject explicitly with "
            "`fasta-lake entrapment inject --stage appended` on the parsimony output."
        )

    edir = Path(entrapment_dir) if entrapment_dir else Path(str(output) + ".entrapment")
    manifest = prepare_from_config(cfg, database, edir)
    click.echo(
        f"[entrapment] class={manifest['class']} stage={manifest['stage']} "
        f"target={manifest['n_target']:,} entrapment={manifest['n_entrapment']:,} "
        f"(ratio {manifest['ratio_entrapment_to_target']})"
    )
    click.echo(f"[entrapment] manifest: {manifest['manifest_path']}")
    return manifest["output_fasta"]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


@cli.command()
@click.argument("config_file", type=click.Path(exists=True))
def validate(config_file: str) -> None:
    """Validate a FASTA Lake JSON configuration file."""
    try:
        cfg = load_config(config_file)
        cfg.validate()
        click.echo("Configuration is valid")
        click.echo(f"  Sources: {len(cfg.database.sources)}")
        click.echo(f"  Strategy: {cfg.inference.strategy}")
        click.echo(f"  Output: {cfg.output.directory}")
    except Exception as e:
        click.echo(f"Configuration error: {e}", err=True)
        raise SystemExit(1)


@cli.command()
@click.option(
    "--template",
    default="full",
    type=click.Choice(["minimal", "clinical", "metaproteomics", "full"]),
    help="Config template type",
)
def init(template: str) -> None:
    """Generate an example fastalake.json configuration."""
    click.echo(create_example_config(template))


# ---------------------------------------------------------------------------
# v5.2: Annotation-aware dedup (Python)
# ---------------------------------------------------------------------------


@cli.command("lake-merge")
@click.option(
    "--source",
    "sources",
    multiple=True,
    required=True,
    metavar="TAG=PATH",
    help="Source FASTA, e.g. --source SP=swissprot.fasta. Repeatable.",
)
@click.option("-o", "--out-fasta", required=True, type=click.Path())
@click.option(
    "-p",
    "--out-provenance",
    default=None,
    type=click.Path(),
    help="Output JSONL sidecar path. Defaults to <out-fasta>.provenance.jsonl",
)
@click.option(
    "--mode",
    "dedup_mode",
    type=click.Choice(["merge", "priority-win"]),
    default="merge",
    help="merge = v5.2 field-level merge (default). "
    "priority-win = legacy behaviour for regression checks.",
)
@click.option(
    "--priority",
    default=None,
    help="Comma-separated tag priority, e.g. UH,SP,TR. Default: SP,TR,UH,GM,GENCODE,GNOMAD,NCBI",
)
@click.option(
    "--embed-sources", is_flag=True, help="Append [FLlake:sources=...] to FASTA headers (opt-in)."
)
@click.option("--compress", is_flag=True, help="Gzip the provenance sidecar.")
@click.option("-v", "--verbose", is_flag=True)
def lake_merge_cmd(
    sources: tuple[str, ...],
    out_fasta: str,
    out_provenance: str | None,
    dedup_mode: str,
    priority: str | None,
    embed_sources: bool,
    compress: bool,
    verbose: bool,
) -> None:
    """Merge multiple source FASTAs into one lake with field-level annotation merge (v5.2).

    Example:

    \b
        fasta-lake lake-merge \\
            --source SP=swissprot.fasta \\
            --source TR=trembl.fasta \\
            --source UH=uhgp.fasta \\
            -o lake.fasta --embed-sources
    """
    _setup_logging(verbose)
    from fasta_lake.lake_builder import build_lake
    from fasta_lake.lake_merge import DEFAULT_PRIORITY

    src_map: dict[str, Path] = {}
    for spec in sources:
        if "=" not in spec:
            raise click.BadParameter(f"--source expects TAG=PATH, got {spec!r}")
        tag, path = spec.split("=", 1)
        src_map[tag.strip()] = Path(path)

    prio = tuple(t.strip() for t in priority.split(",")) if priority else DEFAULT_PRIORITY

    out_fasta_p = Path(out_fasta)
    out_prov_p = (
        Path(out_provenance)
        if out_provenance
        else (
            out_fasta_p.with_suffix(
                out_fasta_p.suffix + (".provenance.jsonl.gz" if compress else ".provenance.jsonl")
            )
        )
    )

    stats = build_lake(
        sources=src_map,
        out_fasta=out_fasta_p,
        out_provenance=out_prov_p,
        dedup_mode=dedup_mode,
        priority=prio,
        embed_sources=embed_sources,
        compress=compress,
    )
    click.echo(f"Lake built: {out_fasta_p}")
    click.echo(f"  Mode:            {stats['dedup_mode']}")
    click.echo(f"  Inputs:          {stats['n_input_sequences']:,}")
    click.echo(
        f"  Unique:          {stats['n_unique_sequences']:,} ({stats['dedup_ratio']}x dedup)"
    )
    click.echo(f"  Shared (≥2 src): {stats['n_shared_across_sources']:,}")
    click.echo(f"  Conflicts:       {stats['n_with_conflicts']:,}")
    click.echo(f"  Ambiguous org:   {stats['n_ambiguous_organism']:,}")
    click.echo(f"  Runtime:         {stats['runtime_seconds']}s")
    click.echo(f"Sidecar: {out_prov_p}")


@cli.command("lake-inspect")
@click.argument("sidecar", type=click.Path(exists=True))
@click.option("--accession", default=None, help="Show record for a specific accession.")
@click.option("--conflicts-only", is_flag=True, help="Only show conflicting records.")
@click.option(
    "--duplicates",
    is_flag=True,
    help="Only show records where ≥2 sources contributed the exact "
    "same sequence (byte-identical collapse).",
)
@click.option("--summary", is_flag=True, help="Print aggregate counts instead of records.")
@click.option("--limit", type=int, default=20, help="Max records to print (default 20).")
def lake_inspect_cmd(
    sidecar: str,
    accession: str | None,
    conflicts_only: bool,
    duplicates: bool,
    summary: bool,
    limit: int,
) -> None:
    """Query a lake provenance sidecar produced by lake-merge."""
    import json as _json
    from collections import Counter

    from fasta_lake.lake_merge import read_provenance_sidecar

    recs = read_provenance_sidecar(Path(sidecar))

    if summary:
        tag_counts: Counter[str] = Counter()
        shared = 0
        conflicts = 0
        ambiguous = 0
        for r in recs:
            for s in r.sources:
                tag_counts[s.tag] += 1
            if len(r.sources) > 1:
                shared += 1
            if r.conflicts:
                conflicts += 1
            if r.organism == "ambiguous":
                ambiguous += 1
        click.echo(f"Total records:     {len(recs):,}")
        click.echo(f"Shared (≥2 src):   {shared:,}")
        click.echo(f"With conflicts:    {conflicts:,}")
        click.echo(f"Ambiguous org:     {ambiguous:,}")
        click.echo("Source contribution:")
        for tag, n in tag_counts.most_common():
            click.echo(f"  {tag:10s} {n:,}")
        return

    shown = 0
    for r in recs:
        if (
            accession
            and r.primary_accession != accession
            and not any(s.accession == accession for s in r.sources)
        ):
            continue
        if conflicts_only and not r.conflicts:
            continue
        if duplicates and len(r.sources) <= 1:
            continue
        click.echo(
            _json.dumps(
                {
                    "primary": f"{r.primary_tag}|{r.primary_accession}",
                    "gene": r.gene,
                    "organism": r.organism,
                    "tax_ids": r.tax_ids,
                    "n_sources": len(r.sources),
                    "alt_accessions": [
                        s.accession for s in r.sources if s.accession != r.primary_accession
                    ],
                    "sources": [f"{s.tag}|{s.accession}" for s in r.sources],
                    "conflicts": r.conflicts,
                },
                indent=2,
            )
        )
        shown += 1
        if shown >= limit:
            click.echo(f"... (limit {limit} reached; use --limit to see more)")
            break


# ---------------------------------------------------------------------------
# FL_MO: multi-omics evidence gate (opt-in; Methods §M5.5)
# ---------------------------------------------------------------------------


@cli.group("fl-mo", hidden=True)
def fl_mo_grp() -> None:
    """Legacy clustered FL_MO reproduction (archived Methods §M5.5).

    Use molecular for the supported exact-sequence TPM workflow.

    Augments per-sample de-novo extraction with metaG/metaT abundance evidence.
    Default (de-novo-only) runs are unaffected; these commands only run when
    explicitly invoked.
    """


@fl_mo_grp.command("gate")
@click.option("--sample", required=True, help="Sample ID")
@click.option(
    "--denovo-fasta",
    required=True,
    type=click.Path(exists=True),
    help="Stage-3 de-novo extracted FASTA (evidence (a))",
)
@click.option(
    "--tpm-tsv",
    required=True,
    type=click.Path(exists=True),
    help="Per-sample TPM TSV with metaG_tpm/metaT_tpm columns",
)
@click.option(
    "--assembly-fasta",
    required=True,
    type=click.Path(exists=True),
    help="Per-sample MEGAHIT predicted-genes FASTA",
)
@click.option(
    "--hash-to-rep",
    required=True,
    type=click.Path(exists=True),
    help="FL_MO hash_to_rep.tsv (sha256 -> cluster rep id)",
)
@click.option(
    "--cluster-rep-fasta",
    required=True,
    type=click.Path(exists=True),
    help="FL_MO cluster97 representative FASTA",
)
@click.option("--output-fasta", required=True, type=click.Path())
@click.option("--output-lookup", required=True, type=click.Path())
@click.option("--output-stats", default=None, type=click.Path())
@click.option(
    "--gate",
    type=click.Choice(["or", "and"]),
    default="or",
    help="or = either evidence rule; and = both rules (no FDR calibration)",
)
@click.option("--metag-top-pct", type=float, default=1.0)
@click.option("--metat-abs", type=float, default=1.0)
def fl_mo_gate(
    sample,
    denovo_fasta,
    tpm_tsv,
    assembly_fasta,
    hash_to_rep,
    cluster_rep_fasta,
    output_fasta,
    output_lookup,
    output_stats,
    gate,
    metag_top_pct,
    metat_abs,
) -> None:
    """Run the FL_MO gate for one sample at Stage 3."""
    from fasta_lake.multi_omics import run_fl_mo_sample

    stats = run_fl_mo_sample(
        sample=sample,
        denovo_fasta=denovo_fasta,
        tpm_tsv=tpm_tsv,
        assembly_fasta=assembly_fasta,
        hash_to_rep=hash_to_rep,
        cluster_rep_fasta=cluster_rep_fasta,
        output_fasta=output_fasta,
        output_lookup=output_lookup,
        output_stats=output_stats,
        gate=gate,
        metag_top_pct=metag_top_pct,
        metat_abs=metat_abs,
    )
    click.echo(
        f"[FL_MO {gate}-gate] sample={stats['sample']} "
        f"denovo={stats['fasta_denovo_written']:,} "
        f"rescued={stats['fasta_rescued_written']:,} "
        f"total={stats['fasta_total']:,} "
        f"(OR={stats['or_gate_members']:,} AND={stats['and_gate_members']:,})"
    )


@fl_mo_grp.command("reproduce")
def fl_mo_reproduce() -> None:
    """Reproduce the published FL_MO razor per-sample median (13,776)."""
    from fasta_lake.multi_omics.reproduce import reproduce

    r = reproduce()
    if not r["matches_published"]:
        raise SystemExit(1)


@cli.command("quantify-study")
@click.option(
    "--groups", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--method",
    type=click.Choice(["sum", "directlfq"]),
    default="sum",
    show_default=True,
    help="Quantify the same fixed study-group assignments.",
)
@click.option("--threads", type=click.IntRange(min=1), default=1, show_default=True)
@click.option(
    "--skip-qc", is_flag=True, help="Explicitly omit automatic QC/PCA for a core-only installation."
)
def quantify_study_cmd(groups, out, method, threads, skip_qc):
    """Quantify completed study groups and retain peptide/ion support diagnostics."""
    import json

    from fasta_lake.quantification import quantify_study

    try:
        if not skip_qc:
            from fasta_lake.downstream import analyze_study, require_analysis

            require_analysis()
        result = quantify_study(groups, out, method=method, threads=threads)
        if not skip_qc:
            result["qc"] = analyze_study(groups, out / "qc", quantification=out)
        else:
            result["qc"] = {"status": "SKIPPED: requested with --skip-qc"}
        (out / "WORKFLOW.json").write_text(json.dumps(result, indent=2) + "\n")
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@cli.command("analyze-study")
@click.option(
    "--groups", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--quantification",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Completed quantify-study output; defaults to the original group sums.",
)
@click.option(
    "--metadata",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="TSV with one unique sample row per acquisition; order is matched by name.",
)
@click.option("--replicate-column", help="Metadata column defining replicate sets; blank excludes.")
@click.option(
    "--pca-samples",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Optional sample roster, one acquisition identifier per line.",
)
@click.option(
    "--pca-features",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Optional fixed group roster; all must be complete within the PCA samples.",
)
def analyze_study_cmd(
    groups, out, quantification, metadata, replicate_column, pca_samples, pca_features
):
    """Export AnnData and run AlphaPeptTools QC, log2 PCA and optional replicate CV."""
    import json

    from fasta_lake.downstream import analyze_study

    try:
        result = analyze_study(
            groups,
            out,
            quantification=quantification,
            metadata=metadata,
            replicate_column=replicate_column,
            pca_samples=pca_samples,
            pca_features=pca_features,
        )
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@cli.command("extract-chunked")
@click.option(
    "--database", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option("--predictions", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--chunk-mib",
    type=click.IntRange(min=1),
    default=256,
    show_default=True,
    help="Reference piece target in MiB; complete records are never split.",
)
@click.option("--chunk-records", type=click.IntRange(min=1), default=500000, show_default=True)
@click.option(
    "--top-fraction",
    type=click.FloatRange(min=0, max=1, min_open=True),
    default=0.30,
    show_default=True,
    help="Rank predictions per sample before matching.",
)
@click.option("--min-length", type=click.IntRange(min=1), default=9, show_default=True)
@click.option("--max-length", type=click.IntRange(min=1), default=50, show_default=True)
@click.option("--normalize-il/--no-normalize-il", default=True, show_default=True)
@click.option("--threads", type=click.IntRange(min=1), default=1, show_default=True)
@click.option("--output-hits", is_flag=True, help="Also retain and merge all peptide hit rows.")
@click.option(
    "--keep-pieces",
    is_flag=True,
    help="Keep selected piece FASTAs/tables after merging; uses additional disk.",
)
def extract_chunked_cmd(
    database,
    predictions,
    out,
    chunk_mib,
    chunk_records,
    top_fraction,
    min_length,
    max_length,
    normalize_il,
    threads,
    output_hits,
    keep_pieces,
):
    """Build evidence in bounded reference pieces and merge before inference."""
    import json

    from fasta_lake.chunked import extract_chunked

    try:
        result = extract_chunked(
            database,
            predictions,
            out,
            chunk_bytes=chunk_mib * 1024 * 1024,
            chunk_records=chunk_records,
            top_fraction=top_fraction,
            min_length=min_length,
            max_length=max_length,
            normalize_il=normalize_il,
            threads=threads,
            output_hits=output_hits,
            keep_pieces=keep_pieces,
        )
    except (OSError, ValueError, RuntimeError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps({"status": result["status"], "stats": result["stats"]}, indent=2))


@cli.command("plan-resources")
@click.option("--memory-budget-gib", type=float)
@click.option("--memory-fraction", type=float, default=0.35, show_default=True)
@click.option("--cpu-fraction", type=float, default=0.5, show_default=True)
@click.option("--threads", type=click.IntRange(min=1))
@click.option("--chunk-mib", type=click.IntRange(min=1))
def plan_resources_cmd(memory_budget_gib, memory_fraction, cpu_fraction, threads, chunk_mib):
    """Show laptop planning choices without starting a workflow."""
    import json

    from fasta_lake.capacity import plan_laptop

    try:
        plan = plan_laptop(
            memory_gib=memory_budget_gib,
            memory_fraction=memory_fraction,
            cpu_fraction=cpu_fraction,
            threads=threads,
            chunk_mib=chunk_mib,
        )
    except (OSError, ValueError, RuntimeError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(plan, indent=2))


from fasta_lake.annotation_cli import annotation  # noqa: E402
from fasta_lake.download_cli import resources  # noqa: E402
from fasta_lake.exploration_cli import explore_cmd  # noqa: E402
from fasta_lake.multi_omics.cli import molecular  # noqa: E402
from fasta_lake.network_cli import network  # noqa: E402
from fasta_lake.taxonomy_cli import taxonomy  # noqa: E402

cli.add_command(molecular)
cli.add_command(annotation)
cli.add_command(resources)
cli.add_command(explore_cmd)
cli.add_command(network)
cli.add_command(taxonomy)


if __name__ == "__main__":
    cli()
