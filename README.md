# NUMTARIS
<img width="1254" height="1254" alt="ChatGPT Image 2026年7月30日 10_57_06" src="https://github.com/user-attachments/assets/388ead86-c204-4e33-98b4-ea80a279cdef" />

# NUMTs Detection Pipeline

Detects candidate NUMTs (nuclear-embedded mitochondrial DNA sequences) from
short-read whole-genome sequencing (WGS) BAM files, using MT-nuclear
discordant- and split-read evidence.

This pipeline builds on the discordant/split-read clustering approach
described in Wei et al., *Nature* 611:105–114 (2022) (see [Citation](#citation)),
reimplemented here as a self-contained, portable pipeline of three scripts.

## Table of contents

- [How it works](#how-it-works)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output format](#output-format)
- [Parameters](#parameters)
- [Reproducibility](#reproducibility)
- [License](#license)

## How it works

```
run_NUMTs_All_fixed.sh   (batch driver, loops over all BAMs in a directory)
        |
        v
NUMTs_detection_fixed.sh  (per-sample: extract MT-related reads, then cluster)
        |
        v
searchNumtCluster_fromDiscordantReads_fixed.py  (cluster + filter -> TSV)
```

1. **`NUMTs_detection_fixed.sh`** — for a single BAM:
   - Indexes the BAM if no `.bai` index exists.
   - Extracts primary-mapped read pairs where either mate maps to the
     mitochondrial contig (`chrM`/`MT`), and splits them into discordant vs.
     split-read SAM files using [samblaster](https://github.com/GregoryFaust/samblaster).
   - Passes those SAM files to `searchNumtCluster_fromDiscordantReads_fixed.py`
     to produce a per-sample TSV of NUMT candidates.

2. **`searchNumtCluster_fromDiscordantReads_fixed.py`** — for a sample's
   discordant/split SAM files:
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

3. **`run_NUMTs_All_fixed.sh`** — runs step 1 for every `*.bam` file in a
   directory, writing each sample's output to its own subdirectory.

> **Scope note:** this pipeline identifies candidate NUMT insertion sites
> from read-pair evidence alone. It does not perform sequence-level
> breakpoint validation (e.g. BLAT realignment of soft-clipped sequence);
> the `filter_status` column in the output should be treated as a
> confidence tier for candidates, not a final call.

## Repository structure

```
.
├── Scripts/
│   ├── NUMTs_detection_fixed.sh                    # per-sample detection script
│   ├── run_NUMTs_All_fixed.sh                      # batch driver over a BAM directory
│   └── searchNumtCluster_fromDiscordantReads_fixed.py  # clustering/filtering logic
├── requirements.txt                            # Python dependencies
├── README.md
└── LICENSE
```

## Requirements

- [samtools](http://www.htslib.org/) (tested with 1.x)
- [samblaster](https://github.com/GregoryFaust/samblaster)
- Python 3.8+ with [pandas](https://pandas.pydata.org/) (see `requirements.txt`)

All three must be available on `PATH`. On macOS with [Homebrew](https://brew.sh/):

```bash
brew install samtools samblaster
```

On Linux with conda:

```bash
conda install -c bioconda samtools samblaster
```

## Installation

```bash
git clone https://github.com/<your-org-or-username>/<repo-name>.git
cd <repo-name>
pip install -r requirements.txt
chmod +x Scripts/NUMTs_detection_fixed.sh Scripts/run_NUMTs_All_fixed.sh
```

No build step is required; the scripts run directly from the repository.

## Test Sample(CCLE Dataset)
URL:https://registry.opendata.aws/depmap-omics-ccle/

## Usage

### Single sample

```bash
./Scripts/NUMTs_detection_fixed.sh <input.bam> <output_dir>
```

Produces, under `<output_dir>`:
- `<sample>.mt.disc.sam` — MT-related discordant read pairs
- `<sample>.mt.split.sam` — MT-related split reads
- `<sample>.NUMTs_candidates.tsv` — clustered NUMT candidates

Example:

```bash
./Scripts/NUMTs_detection_fixed.sh data/sample01.bam results/sample01
```

With non-default thresholds:

```bash
./Scripts/NUMTs_detection_fixed.sh \
    --min-mapq 30 \
    --min-disc-reads 3 \
    --min-split-reads 2 \
    --max-cluster-gap 300 \
    data/sample01.bam results/sample01
```

### Batch (all BAMs in a directory)

```bash
./Scripts/run_NUMTs_All_fixed.sh <bam_dir> <output_root>
```

Runs `NUMTs_detection_fixed.sh` on every `<bam_dir>/*.bam`, writing results to
`<output_root>/<sample>/`.

Example:

```bash
./Scripts/run_NUMTs_All_fixed.sh data/bams/ results/
```

Extra options after the two required arguments (optionally preceded by `--`)
are forwarded to every per-sample run:

```bash
./Scripts/run_NUMTs_All_fixed.sh data/bams/ results/ --min-disc-reads 3
```

### Clustering script directly

`searchNumtCluster_fromDiscordantReads_fixed.py` can also be run standalone
against existing discordant/split SAM files:

```bash
python3 Scripts/searchNumtCluster_fromDiscordantReads_fixed.py \
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

## Output format

`<sample>.NUMTs_candidates.tsv` columns:

| Column                  | Description                                                        |
|--------------------------|----------------------------------------------------------------------|
| `SampleID`              | Sample identifier (BAM filename without extension)                  |
| `chr`                   | Nuclear chromosome of the candidate NUMT                            |
| `start`, `end`          | Candidate NUMT region on `chr` (padded for breakpoint search)       |
| `MT_start`, `MT_end`    | Range of MT coordinates covered by supporting reads                 |
| `NUMT_discordant_reads` | Number of distinct discordant read pairs supporting the cluster     |
| `NUMT_split_reads`      | Number of distinct split reads supporting the cluster               |
| `NUMT_total_reads`      | Sum of discordant + split read support                              |
| `mean_MAPQ`             | Mean mapping quality of the nuclear anchor reads                    |
| `filter_status`         | `PASS_DISC_AND_SPLIT` or `LOW_CONF_DISC_ONLY`                       |
| `discFile`, `splitFile`, `wgsBAM` | Paths to the input files used, for traceability           |

## Parameters

| Parameter          | Flag                  | Default | Meaning                                        |
|----------------------|------------------------|---------|--------------------------------------------------|
| `MIN_MAPQ`          | `--min-mapq`          | 20      | Minimum mapping quality for anchor/split reads    |
| `MIN_DISC_READS`    | `--min-disc-reads`    | 2       | Minimum discordant read pairs to call a candidate |
| `MIN_SPLIT_READS`   | `--min-split-reads`   | 1       | Minimum split reads to mark `PASS_DISC_AND_SPLIT` |
| `MAX_CLUSTER_GAP`   | `--max-cluster-gap`   | 500     | Max bp gap between anchors in the same cluster    |

All four can be overridden from the command line on both
`NUMTs_detection_fixed.sh` and `run_NUMTs_All_fixed.sh` (see [Usage](#usage)).

## Reproducibility

- All detection thresholds are explicit command-line parameters with
  documented defaults (no hidden/hardcoded per-machine paths).
- Script locations are resolved relative to the script's own path
  (`$(dirname "${BASH_SOURCE[0]}")`), so the pipeline runs correctly
  regardless of the working directory or where the repository is cloned.
- Python dependencies are pinned to minimum versions in `requirements.txt`.
- Each output TSV records the exact `discFile`, `splitFile`, and `wgsBAM`
  paths used to generate it, for traceability back to the inputs.

## License

Released under the [MIT License](LICENSE).

