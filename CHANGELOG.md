# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.0.0] - 2026-07-02

### Added
- Initial public release of the NUMT (nuclear-embedded mitochondrial DNA)
  candidate detection pipeline:
  - `NUMTs_detection_fixed.sh` — per-sample extraction of MT-nuclear
    discordant/split reads and clustering into candidate NUMTs.
  - `run_NUMTs_All_fixed.sh` — batch driver over a directory of BAM files.
  - `searchNumtCluster_fromDiscordantReads_fixed.py` — clustering and
    filtering logic producing the candidate TSV.
  - Configurable detection thresholds (`--min-mapq`, `--min-disc-reads`,
    `--min-split-reads`, `--max-cluster-gap`) via command-line flags.
  - README, CITATION, and license for public release.
