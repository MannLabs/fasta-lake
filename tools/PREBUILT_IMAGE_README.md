# FastaLake: ready-to-run Docker demo

This download contains the tested Linux x86_64 image, the complete demo script,
checksums and its exact source commit. Install/open Docker, unpack this download,
then open a terminal in this folder.

```bash
shasum -a 256 -c SHA256SUMS
docker load -i fastalake-docker-linux-amd64.tar.gz
docker tag fastalake:ci fastalake:docker
bash test_container.sh fastalake:docker demo_run_01
```

On Linux, `sha256sum --check SHA256SUMS` is an alternative checksum command.
The script needs Docker and Bash; it does not compile software or download data.
Apple Silicon uses Docker's Linux x86_64 emulation. Each test container has a
4 GiB RAM limit and two CPU slots. Allow disk space for the loaded image and
the new result folder. Use a new output name for each run.

Successful execution ends with `PASS`. Open these files in `demo_run_01/`:

- `campi/study/qc/QC_overview.png`: QC from measured CAMPI acquisition windows.
- `campi/study/study_group_matrix.tsv`: study quantities.
- `campi_directlfq/study_group_matrix.tsv`: directLFQ quantities.
- `synthetic_analysis/sum/qc/PCA.png`: a labelled synthetic PCA teaching example.
- `multiomics/workflow/samples/S01/molecular/gene_evidence.tsv`: molecular selection records.

The [complete demo tour](https://github.com/MannLabs/fasta-lake/blob/main/docs/get-started/demo.md)
explains every example and identifies synthetic inputs. The image includes the
core workflow, directLFQ and AlphaPeptTools. Actual eggNOG databases and optional
ESM-C weights are separate, explicitly chosen resources.

Keep `SOURCE_COMMIT.txt`, `image.json` and the demo logs with your results.
The [first-study tutorial](https://github.com/MannLabs/fasta-lake/blob/main/docs/get-started/first-study.md)
explains how to use your own acquisitions.
