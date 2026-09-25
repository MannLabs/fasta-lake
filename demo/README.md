# FastaLake demo

A self-contained, few-second demonstration of the core idea, runnable on a laptop.

```bash
bash demo/run_demo.sh demo_run_01
```

The script needs Bash, Cargo and a C compiler. It builds the two Rust crates on first use;
the pinned toolchain and dependencies may require network access. Processing takes seconds
after compilation. No Python, spectra or search engine is required. Run through your site's
scheduler if its policy requires this; the demo itself does not require a cluster.

## What it shows

One pooled sequence lake, three samples, and — the point — **three different databases**, each
curated from that sample's own de novo peptide evidence.

```
Stage 0 - build the sequence lake        lake: 20 unique sequences
Stage 2 - evidence filter                evidence lake: 14 of 20 sequences survived
Stage 3 - clustering                     SKIPPED
Stages 4-5 - per-sample curation + parsimony

  SAMPLE       PEPTIDES   DB(seqs)  PARSIMONY
  SAMPLE_A            8         13          4
  SAMPLE_B            8         12          4
  SAMPLE_C            7         10          3
```

Six of the twenty lake sequences carry no supplied prediction peptide and are dropped at stage 2.
Each sample then uses its own synthetic prediction evidence, so `SAMPLE_C` draws 10 sequences
where `SAMPLE_A` has 13 — and they are not subsets of one another. Confirm it yourself:

```bash
diff <(grep '^>' demo_run_01/SAMPLE_A/SAMPLE_A_db.fasta | sort) \
     <(grep '^>' demo_run_01/SAMPLE_C/SAMPLE_C_db.fasta | sort)
```

In a real run the lake holds hundreds of millions of sequences and each sample's database lands
around 5,000–15,000 — but the mechanism is exactly the one above.

The validated workflow uses no sequence clustering. The script retains a legacy MMseqs2
option, but its output is outside the documented expectations above. Leave `MMSEQS` unset
when reproducing this example. Every run requires a fresh output directory.

## The data

Everything here is synthetic and safe to redistribute — roughly 6 KB total.

| File | What it is |
|---|---|
| `lake/reference_lake.fasta` | 20 sequences with headers spanning the formats FastaLake must parse: UniProt `sp\|`, UHGG `MGYG`, GMGC, Prodigal-style assembly, and small-ORF `SMORF`. |
| `predictions/SAMPLE_{A,B,C}_predictions.csv` | De novo peptide predictions, in the format FastaLake consumes: `peptide_prediction_detokenized_unmodified,score`. Overlapping but deliberately distinct subsets. |

No accession is a real database entry, no sequence derives from UHGG, GMGC, GMSC or UniProt bulk
downloads, and nothing here originates from human subject data.

## What it does not cover

This synthetic example ends at database construction. The bundled
[CAMPI example](../examples/campi/README.md) includes measured mzML spectra and SAGE and
continues through search, LFQ, aggregation and study reporting.
