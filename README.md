# NUMTARIS
<img width="1536" height="1024" alt="ツール候補名3" src="https://github.com/user-attachments/assets/00432854-b28e-4452-af51-485a35d4b1ea" />


Inspired by Polaris, the North Star, NUMTARIS aims to serve as a guiding framework for discovering and interpreting mitochondrial DNA insertions in the nuclear genome.
# NUMTs Detection Pipeline

Detects and visualizes candidate NUMTs (nuclear-embedded mitochondrial DNA
segments) from whole-genome sequencing BAM files, by identifying MT-nuclear
discordant and split read pairs and clustering them into candidate
insertion sites.

This pipeline builds on the discordant/split-read clustering approach
described in Wei et al., *Nature* 611:105–114 (2022) (see
[Citation](#citation)), reimplemented here as a self-contained, portable
pipeline.

## Table of contents

- [Pipeline](#pipeline)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Usage](#usage)
- [Output](#output)
- [Parameters](#parameters)
- [Known limitations](#known-limitations)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [License](#license)

## Pipeline

| Script | Role |
|---|---|
| `NUMTs_detection_fixed.sh` | Single-BAM detection: extracts MT-nuclear discordant/split reads (`samtools` + `samblaster`) and clusters them via `numtAnchorCluster.py`. |
| `numtAnchorCluster.py` | Clusters nuclear anchors from the discordant/split SAM files, scores read support, writes `<sample>.NUMTs_candidates.tsv`. |
| `run_NUMTs_All_fixed.sh` | Batch driver: runs detection on every BAM under a directory, then aggregates all samples' TSVs and calls `run_NUMTs_Fig_html.py`. |
| `run_NUMTs_Fig_html.py` | Aggregates `*.NUMTs_candidates.tsv` files (one sample or many) into a circos-style static figure (PNG/SVG/PDF) and an interactive HTML viewer. |

```
run_NUMTs_All_fixed.sh   (batch driver, loops over all BAMs in a directory)
        |
        v
NUMTs_detection_fixed.sh  (per-sample: extract MT-related reads, then cluster)
        |
        v
numtAnchorCluster.py  (cluster + filter -> TSV)
        |
        v
run_NUMTs_Fig_html.py  (aggregate all samples' TSVs -> circos figure + HTML viewer)
```

**`NUMTs_detection_fixed.sh`** — for a single BAM:

1. Indexes the BAM if no `.bai` index exists.
2. Extracts primary-mapped read pairs where either mate maps to the
   mitochondrial contig (`chrM`/`MT`/`M`), and splits them into discordant
   vs. split-read SAM files using `samblaster`.
3. Passes those SAM files to `numtAnchorCluster.py` to produce a per-sample
   TSV of NUMT candidates.

**`numtAnchorCluster.py`** — for a sample's discordant/split SAM files:

- Filters reads by mapping quality (`--min-mapq`).
- For each discordant pair with one mate on the nuclear genome and one on
  MT, takes the nuclear-side coordinate as an anchor.
- Clusters anchors on the same chromosome that are within
  `--max-cluster-gap` bp of each other.
- Keeps clusters supported by at least `--min-disc-reads` distinct read
  pairs, and marks a cluster `PASS_DISC_AND_SPLIT` if it also has at least
  `--min-split-reads` corroborating split reads (otherwise
  `LOW_CONF_DISC_ONLY`).
- Writes one row per candidate NUMT to a TSV.

**`run_NUMTs_All_fixed.sh`** — runs detection for every `*.bam` file in a
directory, writing each sample's output to its own subdirectory, then (unless
`--skip-figures` is passed) calls `run_NUMTs_Fig_html.py` to aggregate all
samples' results.

**`run_NUMTs_Fig_html.py`** — for one or many samples' `*.NUMTs_candidates.tsv`
files:

- Aggregates candidates into two views: `PASS_DISC_AND_SPLIT` only, and
  `PASS_DISC_AND_SPLIT` + `LOW_CONF_DISC_ONLY` combined.
- Renders a circos-style figure (PNG/SVG/PDF) linking mtDNA genes to nuclear
  insertion sites, colored and scaled by read support.
- Renders an interactive Canvas-based HTML viewer with per-sample and
  per-chromosome toggles and hover tooltips.
- Writes ranked TSVs of aggregated links (`NUMTs_top_links.tsv`,
  `NUMTs_top_links_with_mt_gene.tsv`, `NUMTs_raw_records_sorted.tsv`).

**Scope note:** this pipeline identifies candidate NUMT insertion sites from
read-pair evidence alone. It does not perform sequence-level breakpoint
validation (e.g. BLAT realignment of soft-clipped sequence); the
`filter_status` column in the output should be treated as a confidence tier
for candidates, not a final call.

## Repository structure

```
.
├── NUMTs_detection_fixed.sh   # per-sample detection: extract reads + cluster
├── numtAnchorCluster.py       # clustering/filtering logic -> per-sample TSV
├── run_NUMTs_All_fixed.sh     # batch driver over a BAM directory
├── run_NUMTs_Fig_html.py      # aggregation -> circos figure + interactive HTML
├── requirements.txt           # Python dependencies
└── README.md
```

## Requirements

- `samtools` (tested with 1.x)
- `samblaster`
- Python 3.8+ with the packages in `requirements.txt` (`pandas`, `numpy`,
  `matplotlib`):
  ```
  pip install -r requirements.txt
  ```

All three tools must be available on `PATH`. On macOS with Homebrew:

```
brew install samtools samblaster
```

On Linux with conda:

```
conda install -c bioconda samtools samblaster
```

No build step is required beyond making the shell scripts executable:

```
chmod +x NUMTs_detection_fixed.sh run_NUMTs_All_fixed.sh
```

## Usage

### Single sample

```
./NUMTs_detection_fixed.sh [options] input.bam output_dir/
```

Produces, under `output_dir/`:

- `<sample>.mt.disc.sam` — MT-related discordant read pairs
- `<sample>.mt.split.sam` — MT-related split reads
- `<sample>.NUMTs_candidates.tsv` — clustered NUMT candidates

Options: `--min-mapq`, `--min-disc-reads`, `--min-split-reads`,
`--max-cluster-gap` (see `./NUMTs_detection_fixed.sh -h`).

Example with non-default thresholds:

```
./NUMTs_detection_fixed.sh \
    --min-mapq 30 \
    --min-disc-reads 3 \
    --min-split-reads 2 \
    --max-cluster-gap 300 \
    data/sample01.bam results/sample01
```

### Batch (all BAMs in a directory)

```
./run_NUMTs_All_fixed.sh bam_dir/ output_root/ [--skip-figures] [-- extra detection options]
```

Runs detection on every `*.bam` in `bam_dir/`, writing each sample's output
to `output_root/<sample>/`. Unless `--skip-figures` is passed, it then
aggregates every sample's TSV under `output_root/` and writes figures to
`output_root/output_<label>/`.

Example overriding detection thresholds for every sample:

```
./run_NUMTs_All_fixed.sh bam_dir/ output_root/ -- --min-disc-reads 3
```

### Figures only

Useful for re-rendering after the fact, or for aggregating TSVs that were
generated separately:

```
python3 run_NUMTs_Fig_html.py --input-dir output_root/ [--output-dir output_root/]
```

`--input-dir` is searched recursively for `*.NUMTs_candidates.tsv`, so it
accepts either a single sample's output directory or a root directory
containing one subdirectory per sample.

### Clustering script directly

`numtAnchorCluster.py` can also be run standalone against existing
discordant/split SAM files:

```
python3 numtAnchorCluster.py \
    --sample sample01 \
    --bam data/sample01.bam \
    --disc results/sample01/sample01.mt.disc.sam \
    --split results/sample01/sample01.mt.split.sam \
    --out results/sample01/sample01.NUMTs_candidates.tsv \
    --min-mapq 20 \
    --min-disc-reads 2 \
    --min-split-reads 1 \
    --max-cluster-gap 500
```

## Output

### Per-sample TSV

`<sample>.NUMTs_candidates.tsv` columns:

| Column | Description |
|---|---|
| `SampleID` | Sample identifier (BAM filename without extension) |
| `chr` | Nuclear chromosome of the candidate NUMT |
| `start`, `end` | Candidate NUMT region on `chr` (padded for breakpoint search) |
| `MT_start`, `MT_end` | Range of MT coordinates covered by supporting reads |
| `NUMT_discordant_reads` | Number of distinct discordant read pairs supporting the cluster |
| `NUMT_split_reads` | Number of distinct split reads supporting the cluster |
| `NUMT_total_reads` | Deduplicated count of discordant + split supporting read pairs |
| `mean_MAPQ` | Mean mapping quality of the nuclear anchor reads |
| `filter_status` | `PASS_DISC_AND_SPLIT` or `LOW_CONF_DISC_ONLY` |
| `discFile`, `splitFile`, `wgsBAM` | Paths to the input files used, for traceability |

### Aggregated figures and tables

`run_NUMTs_Fig_html.py` aggregates candidates across samples into two views
(`PASS_DISC_AND_SPLIT` only, and `PASS_DISC_AND_SPLIT` +
`LOW_CONF_DISC_ONLY`), each written to its own `output_<label>/` folder
containing:

- `NUMTs_semicircle_rotated90_segmented_chrM_grayRing_coloredLinks.{png,svg,pdf}` — circos-style static figure linking mtDNA genes to nuclear insertion sites
- `NUMTs_interactive_viewer.html` — interactive Canvas-based viewer with per-sample/per-chromosome toggles and hover tooltips
- `NUMTs_top_links.tsv` — aggregated links ranked by total supporting reads
- `NUMTs_top_links_with_mt_gene.tsv` — links grouped by mtDNA gene and nuclear chromosome
- `NUMTs_raw_records_sorted.tsv` — all contributing per-sample records, sorted by read support

## Parameters

| Parameter | Flag | Default | Meaning |
|---|---|---|---|
| `MIN_MAPQ` | `--min-mapq` | 20 | Minimum mapping quality for anchor/split reads |
| `MIN_DISC_READS` | `--min-disc-reads` | 2 | Minimum discordant read pairs to call a candidate |
| `MIN_SPLIT_READS` | `--min-split-reads` | 1 | Minimum split reads to mark `PASS_DISC_AND_SPLIT` |
| `MAX_CLUSTER_GAP` | `--max-cluster-gap` | 500 | Max bp gap between anchors in the same cluster |

All four can be overridden from the command line on both
`NUMTs_detection_fixed.sh` and `run_NUMTs_All_fixed.sh` (forwarded via
`-- <options>` in the latter).

## Known limitations

- Chromosome-naming logic (`is_primary_chrom` in `numtAnchorCluster.py`;
  `chr_len` / `normalize_chr` in `run_NUMTs_Fig_html.py`) assumes
  UCSC/Ensembl-style GRCh38 contig names (`chr1`..`chr22`, `chrX`, `chrY`,
  `chrM`/`MT`). RefSeq accession-style names (e.g. `NC_000001.11`) are not
  supported and will produce empty results.
- `NUMTs_detection_fixed.sh` uses `samtools view -F 2` to drop
  "properly paired" alignment records before running `samblaster`. On some
  aligners this flag can also be set on supplementary (split) alignment
  records belonging to an otherwise well-behaved read pair, which could
  suppress genuine split-read evidence. If split-read counts look
  unexpectedly low on your data, this is worth validating against known
  positive NUMTs.

## Reproducibility

- All detection thresholds are explicit command-line parameters with
  documented defaults (no hidden/hardcoded per-machine paths).
- Script locations are resolved relative to the script's own path
  (`$(dirname "${BASH_SOURCE[0]}")`), so the pipeline runs correctly
  regardless of the working directory or where the repository is cloned.
- Python dependencies are pinned to minimum versions in `requirements.txt`.
- Each output TSV records the exact `discFile`, `splitFile`, and `wgsBAM`
  paths used to generate it, for traceability back to the inputs.

## Citation

This pipeline reimplements the discordant/split-read clustering approach
for NUMT detection described in:

> Wei, W., Schon, K.R., Elgar, G. et al. Nuclear-embedded mitochondrial DNA
> sequences in 66,083 human genomes. *Nature* 611, 105–114 (2022).
> https://doi.org/10.1038/s41586-022-05288-7

## License

Released under the MIT License.

