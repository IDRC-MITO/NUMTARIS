# Contributing

Thanks for your interest in improving this NUMT detection pipeline.

## Reporting issues

Please open a GitHub issue and include:
- The command you ran (with any non-default thresholds).
- The version/commit of the pipeline you used.
- `samtools --version` and `samblaster --version` output.
- Relevant log output (the pipeline's own `echo`/logging messages are usually enough — please do not attach BAM/SAM files or any sequencing data).

## Submitting changes

1. Fork the repository and create a branch for your change.
2. Keep pull requests focused on a single change; avoid bundling unrelated fixes.
3. For shell scripts: keep `set -euo pipefail`, quote variable expansions, and prefer explicit `--long-flag` options over positional magic.
4. For Python: follow [PEP 8](https://peps.python.org/pep-0008/) and keep functions small and testable.
5. Update `README.md` and `CHANGELOG.md` if your change affects usage, output format, or default parameters.
6. Do not commit any real sequencing data, sample identifiers, or file paths from your own environment — only synthetic/example data belongs in the repository.

## Development notes

There is no formal test suite yet. At minimum, before submitting a change:
- Run `bash -n <script>.sh` on any modified shell script.
- Run `python3 -m py_compile <script>.py` on any modified Python script.
- Run the modified script against a small test BAM (e.g. a few thousand reads) end to end and confirm the output TSV has the expected columns.
