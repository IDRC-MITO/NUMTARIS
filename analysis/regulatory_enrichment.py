#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regulatory-element enrichment and TSS-proximity analysis for NUMT calls.

Given the per-sample NUMT call tables produced by the NUMTARIS detection
pipeline (``*.NUMTs_candidates.tsv``), this script computes, per cohort:

  1. Enrichment of high-confidence NUMTs at each supplied regulatory
     annotation (e.g. Fanta.Bio CREs, FANTOM5 enhancers, FANTOM5
     CAGE-defined promoters, ENCODE cCREs) relative to chromosome-,
     length- and mappability-matched random genomic intervals.
  2. For one designated "promoter" annotation (default FANTOM5 CAGE
     peaks): the fold enrichment across cohorts and the observed vs.
     expected number of overlapping NUMTs.
  3. The distance from each NUMT interval midpoint to the nearest
     protein-coding gene TSS (GENCODE), observed vs. matched random
     intervals.

Method
------
* Only ``PASS_DISC_AND_SPLIT`` (high-confidence) NUMT calls are used.
* NUMT calls are pooled across all samples within a cohort and are not
  de-duplicated (a recurrent locus contributes once per carrier).
* For each NUMT, matched random intervals are drawn on the *same*
  chromosome, preserving the observed interval length, and (when a Umap
  hg38 k=100 mappability track is supplied) restricted to fully mappable
  regions (mappability score == 1.0).
* ``--n-random`` (default 1000) matched random datasets are generated
  with a fixed seed (``--seed``, default 12345).
* Enrichment = observed overlap fraction / mean overlap fraction across
  random datasets.
* Empirical P = (k + 1) / (N + 1), where k is the number of random
  datasets with a value at least as extreme as the observed value and
  N = ``--n-random``.

Inputs
------
--numt-root DIR   Directory containing one sub-directory per cohort.
                  Each cohort sub-directory holds that cohort's
                  ``*.NUMTs_candidates.tsv`` files (one per sample).
--ref-dir DIR     Directory containing the annotation BED files, the
                  GENCODE GTF and (optionally) the Umap bedGraph.
--config FILE     Optional JSON overriding cohort labels / order and
                  annotation file names (see config.example.json).
--out-dir DIR     Output directory (created if missing).

External reference files (place under --ref-dir, names set in --config)
----------------------------------------------------------------------
  <regulatory annotation BEDs>   e.g. FANTOM5, ENCODE SCREEN cCRE, Fanta.Bio
  gencode.v49.annotation.gtf[.gz]   GENCODE  (https://www.gencodegenes.org)
  hg38.umap.k100.bedGraph          Umap k=100 mappability (https://bismap.hoffmanlab.org) [optional]
  hg38.chrom.sizes                 UCSC hg38 chrom sizes [optional; built-in fallback]

Outputs
-------
  regulatory_enrichment_summary.pdf / .png    combined multi-panel figure
  enrichment_by_annotation.{pdf,png}
  promoter_enrichment_by_cohort.{pdf,png}
  promoter_observed_vs_expected.{pdf,png}
  nearest_tss_distance.{pdf,png}
  source_data/annotation_enrichment.tsv
  source_data/promoter_enrichment_by_cohort.tsv
  source_data/promoter_observed_vs_expected.tsv
  source_data/nearest_tss_distance.tsv
  source_data/numts_used.tsv
  source_data/random_datasets_annotation_overlap_counts.tsv
  source_data/random_datasets_nearest_tss_per_dataset.tsv

Dependencies: python>=3.9, numpy, pandas, matplotlib.
"""

import argparse
import bisect
import csv
import glob
import gzip
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Defaults (override with --config)
# ----------------------------------------------------------------------
FILTER_STATUS = "PASS_DISC_AND_SPLIT"

DEFAULT_CONFIG = {
    # cohort sub-directory (under --numt-root)  ->  display label
    "cohorts": {
        "MS00": "Control",
        "MS01": "Psychiatry",
        "MS03": "Esophageal cancer",
        "MS04": "Pancreatic cancer",
        "MS05": "Colorectal cancer",
        "Mito_Blood": "Mitochondrial disease",
    },
    # left-to-right plotting order
    "cohort_order": [
        "Control",
        "Psychiatry",
        "Esophageal cancer",
        "Pancreatic cancer",
        "Colorectal cancer",
        "Mitochondrial disease",
    ],
    "cohort_colors": {
        "Control": "#9e9e9e",
        "Psychiatry": "#7fc7e8",
        "Esophageal cancer": "#2e7d32",
        "Pancreatic cancer": "#00695c",
        "Colorectal cancer": "#e8622c",
        "Mitochondrial disease": "#d81b8c",
    },
    # annotation key -> BED file name under --ref-dir
    "annotation_files": {
        "Fanta.Bio CRE": "human-CREv1.0.hg38.cre-peaks.CREname.bed",
        "FANTOM5 enhancer": "F5.hg38.enhancers.bed",
        "FANTOM5 CAGE peak": "hg38_fair+new_CAGE_peaks_phase1and2.bed",
        "ENCODE PLS": "GRCh38-PLS.bed",
        "ENCODE ELS": "GRCh38-ELS.bed",
        "ENCODE CTCF": "GRCh38-CTCF.bed",
    },
    # left-to-right order of annotation groups
    "annotation_order": [
        "ENCODE CTCF",
        "ENCODE ELS",
        "ENCODE PLS",
        "FANTOM5 CAGE peak",
        "FANTOM5 enhancer",
        "Fanta.Bio CRE",
    ],
    # the annotation used for the per-cohort and observed-vs-expected plots
    "promoter_annotation": "FANTOM5 CAGE peak",
    "gencode_gtf": "gencode.v49.annotation.gtf",
    "umap_bedgraph": "hg38.umap.k100.bedGraph",
    "chrom_sizes": "hg38.chrom.sizes",
}

# GRCh38 primary assembly chromosome sizes (used if hg38.chrom.sizes is absent)
DEFAULT_HG38_CHROM_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------
def log(msg):
    print(msg, flush=True)


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def safe_int(x):
    return int(float(str(x).replace(",", "").strip()))


def detect_delimiter(path):
    return "," if str(path).endswith(".csv") else "\t"


def significance_stars(p):
    try:
        p = float(p)
    except (TypeError, ValueError):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


# ----------------------------------------------------------------------
# Interval set with fast overlap queries (merged, sorted)
# ----------------------------------------------------------------------
class IntervalSet:
    """Merged, sorted non-overlapping intervals per chromosome."""

    def __init__(self):
        self._starts = {}
        self._ends = {}

    @classmethod
    def from_pairs(cls, pairs_by_chr):
        obj = cls()
        for chrom, pairs in pairs_by_chr.items():
            merged = _merge_intervals(pairs)
            if not merged:
                continue
            obj._starts[chrom] = np.array([s for s, _ in merged], dtype=np.int64)
            obj._ends[chrom] = np.array([e for _, e in merged], dtype=np.int64)
        return obj

    def chromosomes(self):
        return set(self._starts)

    def overlaps(self, chrom, start, end):
        starts = self._starts.get(chrom)
        if starts is None:
            return False
        ends = self._ends[chrom]
        # rightmost interval whose start < end
        j = bisect.bisect_right(starts, end - 1) - 1
        if j >= 0 and ends[j] > start:
            return True
        # an interval starting within [start, end)
        k = bisect.bisect_left(starts, start)
        if k < len(starts) and starts[k] < end:
            return True
        return False


def _merge_intervals(pairs):
    pairs = sorted((int(s), int(e)) for s, e in pairs if e > s)
    if not pairs:
        return []
    merged = [list(pairs[0])]
    for s, e in pairs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


# ----------------------------------------------------------------------
# Mappability-matched random sampler
# ----------------------------------------------------------------------
class RandomIntervalSampler:
    """
    Draw a random interval of a given length on a given chromosome.

    If a mappable-interval set is supplied, the interval is drawn so that
    its start falls on a mappable base and it is accepted only if it fits
    entirely inside one mappable interval.  Otherwise the interval is
    drawn uniformly across the chromosome (chromosome- and length-matched
    only).
    """

    def __init__(self, chrom_sizes, mappable_pairs_by_chr, rng):
        self.chrom_sizes = chrom_sizes
        self.rng = rng
        self._starts = {}
        self._ends = {}
        self._cum = {}  # cumulative mappable length, prefix-exclusive
        self._total = {}
        if mappable_pairs_by_chr:
            for chrom, pairs in mappable_pairs_by_chr.items():
                merged = _merge_intervals(pairs)
                if not merged:
                    continue
                starts = np.array([s for s, _ in merged], dtype=np.int64)
                ends = np.array([e for _, e in merged], dtype=np.int64)
                lengths = ends - starts
                cum = np.concatenate([[0], np.cumsum(lengths)])
                self._starts[chrom] = starts
                self._ends[chrom] = ends
                self._cum[chrom] = cum
                self._total[chrom] = int(cum[-1])

    def has_mappable(self, chrom):
        return chrom in self._total and self._total[chrom] > 0

    def draw(self, chrom, length, max_try=100):
        length = int(length)
        chrom_len = self.chrom_sizes.get(chrom)
        if chrom_len is None or chrom_len <= length:
            return None

        if self.has_mappable(chrom):
            starts = self._starts[chrom]
            ends = self._ends[chrom]
            cum = self._cum[chrom]
            total = self._total[chrom]
            for _ in range(max_try):
                k = self.rng.randrange(total)
                idx = bisect.bisect_right(cum, k) - 1
                pos = int(starts[idx] + (k - cum[idx]))
                if ends[idx] - pos >= length:
                    return chrom, pos, pos + length
            return None  # could not place within a mappable interval

        start = self.rng.randint(0, chrom_len - length)
        return chrom, start, start + length


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------
def load_chrom_sizes(path):
    if path and Path(path).exists():
        sizes = {}
        with open(path) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2 and cols[0].strip():
                    try:
                        sizes[cols[0]] = int(cols[1])
                    except ValueError:
                        continue
        if sizes:
            log(f"[LOAD] chrom sizes: {path} ({len(sizes)} contigs)")
            return sizes
    log("[INFO] using built-in GRCh38 primary chromosome sizes")
    return dict(DEFAULT_HG38_CHROM_SIZES)


def load_bed_pairs(path, valid_chroms):
    pairs = defaultdict(list)
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line[0] == "#" or line.startswith(("track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3 or cols[0] not in valid_chroms:
                continue
            try:
                s, e = safe_int(cols[1]), safe_int(cols[2])
            except ValueError:
                continue
            if e > s:
                pairs[cols[0]].append((s, e))
    return pairs


def load_umap_pairs(path, valid_chroms, min_score=1.0):
    if not path or not Path(path).exists():
        log(f"[WARN] Umap track not found ({path}); random intervals will be "
            f"chromosome- and length-matched only")
        return None
    pairs = defaultdict(list)
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line[0] == "#" or line.startswith(("track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4 or cols[0] not in valid_chroms:
                continue
            try:
                s, e, score = safe_int(cols[1]), safe_int(cols[2]), float(cols[3])
            except ValueError:
                continue
            if e > s and score >= min_score:
                pairs[cols[0]].append((s, e))
    n = sum(len(v) for v in pairs.values())
    log(f"[LOAD] Umap mappable intervals (score >= {min_score}): {n}")
    return pairs


def load_protein_coding_tss(gtf_path, valid_chroms):
    if not Path(gtf_path).exists():
        raise FileNotFoundError(f"GENCODE GTF not found: {gtf_path}")
    log(f"[LOAD] GENCODE GTF: {gtf_path}")
    tss_by_chr = defaultdict(list)
    with open_text(gtf_path) as f:
        for line in f:
            if not line.strip() or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene" or cols[0] not in valid_chroms:
                continue
            attrs = {}
            for item in cols[8].strip().split(";"):
                item = item.strip()
                if " " in item:
                    key, val = item.split(" ", 1)
                    attrs[key] = val.replace('"', "").strip()
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))
            if gene_type != "protein_coding":
                continue
            strand = cols[6]
            try:
                start, end = safe_int(cols[3]), safe_int(cols[4])
            except ValueError:
                continue
            if strand == "+":
                tss = start - 1
            elif strand == "-":
                tss = end - 1
            else:
                continue
            tss_by_chr[cols[0]].append(tss)
    tss_arrays = {c: np.array(sorted(v), dtype=np.int64) for c, v in tss_by_chr.items()}
    log(f"[OK] protein-coding TSS: {sum(len(v) for v in tss_arrays.values())}")
    return tss_arrays


def nearest_tss_distance(chrom, start, end, tss_arrays):
    arr = tss_arrays.get(chrom)
    if arr is None or len(arr) == 0:
        return None
    midpoint = (start + end) // 2
    i = bisect.bisect_left(arr, midpoint)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(arr):
            d = abs(int(arr[j]) - midpoint)
            if best is None or d < best:
                best = d
    return best


# ----------------------------------------------------------------------
# Read NUMT calls
# ----------------------------------------------------------------------
def read_numt_calls(numt_root, cohorts, valid_chroms):
    rows = []
    for subdir, label in cohorts.items():
        cohort_dir = Path(numt_root) / subdir
        if not cohort_dir.is_dir():
            log(f"[SKIP] missing cohort directory: {cohort_dir}")
            continue
        files = []
        for ext in ("*.NUMTs_candidates.tsv", "*.tsv", "*.txt", "*.csv"):
            files.extend(sorted(glob.glob(str(cohort_dir / ext))))
        files = sorted(set(files))
        n_before = len(rows)
        for path in files:
            delim = detect_delimiter(path)
            with open(path, newline="") as f:
                reader = csv.DictReader(f, delimiter=delim)
                needed = {"SampleID", "chr", "start", "end", "filter_status"}
                if not needed.issubset(set(reader.fieldnames or [])):
                    continue
                for r in reader:
                    if r.get("filter_status") != FILTER_STATUS:
                        continue
                    chrom = r["chr"]
                    if chrom not in valid_chroms:
                        continue
                    try:
                        s, e = safe_int(r["start"]), safe_int(r["end"])
                    except ValueError:
                        continue
                    if e <= s:
                        continue
                    rows.append({
                        "cohort_dir": subdir,
                        "Cohort": label,
                        "SampleID": r.get("SampleID", "."),
                        "chr": chrom,
                        "start": s,
                        "end": e,
                        "length": e - s,
                        "source_file": Path(path).name,
                    })
        log(f"[OK] {label} ({subdir}): {len(rows) - n_before} high-confidence NUMTs")
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Core: one randomization pass serving all outputs
# ----------------------------------------------------------------------
def run_analysis(numts, annotations, tss_arrays, sampler, n_random, cohort_order):
    ann_keys = list(annotations)
    numt_records = numts.to_dict("records")

    # ---- observed ----
    obs_overlap = {c: {a: 0 for a in ann_keys} for c in cohort_order}
    obs_n = {c: 0 for c in cohort_order}
    obs_tss = defaultdict(list)

    for r in numt_records:
        c = r["Cohort"]
        if c not in obs_n:
            continue
        obs_n[c] += 1
        for a in ann_keys:
            if annotations[a].overlaps(r["chr"], r["start"], r["end"]):
                obs_overlap[c][a] += 1
        d = nearest_tss_distance(r["chr"], r["start"], r["end"], tss_arrays)
        if d is not None:
            obs_tss[c].append(d)

    # ---- random datasets ----
    rand_overlap_counts = {c: {a: [] for a in ann_keys} for c in cohort_order}
    rand_tss_median = {c: [] for c in cohort_order}
    rand_tss_mean = {c: [] for c in cohort_order}

    for it in range(n_random):
        if (it + 1) % max(1, n_random // 10) == 0:
            log(f"[RANDOM] dataset {it + 1}/{n_random}")
        set_counts = {c: {a: 0 for a in ann_keys} for c in cohort_order}
        set_tss = defaultdict(list)
        for r in numt_records:
            c = r["Cohort"]
            if c not in obs_n:
                continue
            drawn = sampler.draw(r["chr"], r["length"])
            if drawn is None:
                continue
            rc, rs, re = drawn
            for a in ann_keys:
                if annotations[a].overlaps(rc, rs, re):
                    set_counts[c][a] += 1
            d = nearest_tss_distance(rc, rs, re, tss_arrays)
            if d is not None:
                set_tss[c].append(d)
        for c in cohort_order:
            for a in ann_keys:
                rand_overlap_counts[c][a].append(set_counts[c][a])
            if set_tss[c]:
                rand_tss_median[c].append(float(np.median(set_tss[c])))
                rand_tss_mean[c].append(float(np.mean(set_tss[c])))

    # ---- enrichment table ----
    enrich_rows = []
    for c in cohort_order:
        n = obs_n[c]
        if n == 0:
            continue
        for a in ann_keys:
            obs_hits = obs_overlap[c][a]
            rc = np.array(rand_overlap_counts[c][a], dtype=float)
            rand_mean = float(rc.mean()) if rc.size else float("nan")
            rand_sd = float(rc.std(ddof=0)) if rc.size else float("nan")
            obs_frac = obs_hits / n
            rand_frac = rand_mean / n if n else float("nan")
            fold = obs_frac / rand_frac if rand_frac and rand_frac > 0 else float("nan")
            k = int(np.sum(rc >= obs_hits)) if rc.size else 0
            emp_p = (k + 1) / (rc.size + 1) if rc.size else float("nan")
            enrich_rows.append({
                "Cohort": c,
                "Annotation": a,
                "N_NUMTs": n,
                "Observed_overlap": obs_hits,
                "Observed_overlap_fraction": obs_frac,
                "Expected_overlap_mean": rand_mean,
                "Expected_overlap_sd": rand_sd,
                "Enrichment_fold": fold,
                "Empirical_P": emp_p,
                "Significance": significance_stars(emp_p),
                "N_random": int(rc.size),
            })
    enrich_df = pd.DataFrame(enrich_rows)

    # ---- nearest-TSS table ----
    tss_rows = []
    for c in cohort_order:
        od = np.array(obs_tss[c], dtype=float)
        rm = np.array(rand_tss_median[c], dtype=float)
        rmean = np.array(rand_tss_mean[c], dtype=float)
        if od.size == 0 or rm.size == 0:
            continue
        obs_median = float(np.median(od))
        obs_mean = float(np.mean(od))
        k_med = int(np.sum(rm <= obs_median))
        k_mean = int(np.sum(rmean <= obs_mean))
        tss_rows.append({
            "Cohort": c,
            "N_NUMTs": int(od.size),
            "Observed_median_TSS_distance_bp": obs_median,
            "Expected_median_TSS_distance_bp_mean": float(np.mean(rm)),
            "Expected_median_TSS_distance_bp_sd": float(np.std(rm, ddof=0)),
            "Observed_mean_TSS_distance_bp": obs_mean,
            "Expected_mean_TSS_distance_bp_mean": float(np.mean(rmean)),
            "Expected_mean_TSS_distance_bp_sd": float(np.std(rmean, ddof=0)),
            "Empirical_P_median_observed_closer": (k_med + 1) / (rm.size + 1),
            "Empirical_P_mean_observed_closer": (k_mean + 1) / (rmean.size + 1),
            "Significance": significance_stars((k_med + 1) / (rm.size + 1)),
            "N_random": int(rm.size),
        })
    tss_df = pd.DataFrame(tss_rows)

    # ---- raw per-dataset tables ----
    rand_overlap_long = []
    for c in cohort_order:
        for a in ann_keys:
            for i, v in enumerate(rand_overlap_counts[c][a]):
                rand_overlap_long.append(
                    {"Cohort": c, "Annotation": a, "Random_dataset": i, "Overlap_count": v}
                )
    rand_overlap_long_df = pd.DataFrame(rand_overlap_long)

    rand_tss_long = []
    for c in cohort_order:
        for i, (md, mn) in enumerate(zip(rand_tss_median[c], rand_tss_mean[c])):
            rand_tss_long.append(
                {"Cohort": c, "Random_dataset": i,
                 "Median_TSS_distance_bp": md, "Mean_TSS_distance_bp": mn}
            )
    rand_tss_long_df = pd.DataFrame(rand_tss_long)

    return enrich_df, tss_df, rand_overlap_long_df, rand_tss_long_df


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
def plot_enrichment_by_annotation(ax, enrich_df, cfg):
    order_ann = [a for a in cfg["annotation_order"] if a in set(enrich_df["Annotation"])]
    order_coh = [c for c in cfg["cohort_order"] if c in set(enrich_df["Cohort"])]
    pivot = (enrich_df.pivot(index="Annotation", columns="Cohort", values="Enrichment_fold")
             .reindex(index=order_ann, columns=order_coh))
    x = np.arange(len(order_ann))
    w = 0.8 / max(1, len(order_coh))
    for i, c in enumerate(order_coh):
        ax.bar(x + (i - (len(order_coh) - 1) / 2) * w, pivot[c].values, w,
               label=c, color=cfg["cohort_colors"].get(c, None))
    ax.axhline(1.0, ls="--", lw=1, color="#1565c0")
    ax.set_xticks(x)
    ax.set_xticklabels(order_ann, rotation=25, ha="right")
    ax.set_ylabel("Observed / random enrichment (fold)")
    ax.set_title("Regulatory-element enrichment of high-confidence NUMTs", loc="left")
    ax.legend(frameon=False, fontsize=7, ncol=2)


def plot_promoter_enrichment_by_cohort(ax, enrich_df, cfg):
    key = cfg["promoter_annotation"]
    sub = enrich_df[enrich_df["Annotation"] == key].copy().sort_values("Enrichment_fold")
    y = np.arange(len(sub))
    ax.barh(y, sub["Enrichment_fold"].values,
            color=[cfg["cohort_colors"].get(c, "#666666") for c in sub["Cohort"]])
    ax.axvline(1.0, ls="--", lw=1, color="#555555")
    ax.set_yticks(y)
    ax.set_yticklabels(sub["Cohort"])
    ax.set_xlabel("Observed / random enrichment (fold)")
    ax.set_title(f"{key} enrichment by cohort", loc="left")
    xmax = float(np.nanmax(sub["Enrichment_fold"].values)) if len(sub) else 1.0
    for i, (_, r) in enumerate(sub.iterrows()):
        ax.text(r["Enrichment_fold"] + xmax * 0.02, i, r["Significance"],
                va="center", ha="left", fontsize=10)
    ax.set_xlim(0, xmax * 1.15)


def plot_promoter_observed_vs_expected(ax, enrich_df, cfg):
    key = cfg["promoter_annotation"]
    sub = enrich_df[enrich_df["Annotation"] == key].copy()
    order_coh = [c for c in cfg["cohort_order"] if c in set(sub["Cohort"])]
    sub = sub.set_index("Cohort").reindex(order_coh)
    x = np.arange(len(order_coh))
    w = 0.38
    ax.bar(x - w / 2, sub["Observed_overlap"].values, w, label="Observed",
           color="#1f77b4")
    ax.bar(x + w / 2, sub["Expected_overlap_mean"].values, w,
           yerr=sub["Expected_overlap_sd"].values, capsize=3,
           label="Matched random (mean)", color="#ff7f0e")
    ax.set_xticks(x)
    ax.set_xticklabels(order_coh, rotation=25, ha="right")
    ax.set_ylabel(f"NUMTs overlapping\n{key} (n)")
    ax.set_title("Observed vs. expected promoter-annotation overlap", loc="left")
    ax.legend(frameon=False, fontsize=7)


def plot_nearest_tss_distance(ax, tss_df, cfg):
    order_coh = [c for c in cfg["cohort_order"] if c in set(tss_df["Cohort"])]
    df = tss_df.set_index("Cohort").reindex(order_coh)
    x = np.arange(len(order_coh))
    w = 0.38
    ax.bar(x - w / 2, df["Observed_median_TSS_distance_bp"].values / 1000, w,
           label="Observed", color="#1f77b4")
    ax.bar(x + w / 2, df["Expected_median_TSS_distance_bp_mean"].values / 1000, w,
           yerr=df["Expected_median_TSS_distance_bp_sd"].values / 1000, capsize=3,
           label="Matched random", color="#ff7f0e")
    ax.set_xticks(x)
    ax.set_xticklabels(order_coh, rotation=25, ha="right")
    ax.set_ylabel("Median nearest TSS distance (kb)")
    ax.set_title("Distance to nearest protein-coding TSS", loc="left")
    ax.legend(frameon=False, fontsize=7)
    ymax = np.nanmax([
        np.nanmax(df["Observed_median_TSS_distance_bp"].values / 1000),
        np.nanmax(df["Expected_median_TSS_distance_bp_mean"].values / 1000),
    ])
    for i, (_, r) in enumerate(df.iterrows()):
        ax.text(i, ymax * 1.03, r["Significance"], ha="center", va="bottom", fontsize=10)


def make_figures(enrich_df, tss_df, cfg, out_dir):
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    plot_enrichment_by_annotation(axes[0, 0], enrich_df, cfg)
    plot_promoter_enrichment_by_cohort(axes[0, 1], enrich_df, cfg)
    plot_promoter_observed_vs_expected(axes[1, 0], enrich_df, cfg)
    plot_nearest_tss_distance(axes[1, 1], tss_df, cfg)
    fig.suptitle("NUMT enrichment at regulatory elements and proximity to TSS",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / "regulatory_enrichment_summary.pdf")
    fig.savefig(out_dir / "regulatory_enrichment_summary.png", dpi=300)
    plt.close(fig)

    for name, fn in [
        ("enrichment_by_annotation",
         lambda ax: plot_enrichment_by_annotation(ax, enrich_df, cfg)),
        ("promoter_enrichment_by_cohort",
         lambda ax: plot_promoter_enrichment_by_cohort(ax, enrich_df, cfg)),
        ("promoter_observed_vs_expected",
         lambda ax: plot_promoter_observed_vs_expected(ax, enrich_df, cfg)),
        ("nearest_tss_distance",
         lambda ax: plot_nearest_tss_distance(ax, tss_df, cfg)),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4.6))
        fn(ax)
        fig.tight_layout()
        fig.savefig(out_dir / f"{name}.pdf")
        fig.savefig(out_dir / f"{name}.png", dpi=300)
        plt.close(fig)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Regulatory-element enrichment and TSS-proximity analysis "
                    "for NUMTARIS NUMT calls.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--numt-root", required=True,
                   help="directory with one sub-directory per cohort of "
                        "*.NUMTs_candidates.tsv files")
    p.add_argument("--ref-dir", required=True,
                   help="directory with annotation BEDs, GENCODE GTF and Umap track")
    p.add_argument("--out-dir", required=True, help="output directory")
    p.add_argument("--config", default=None,
                   help="optional JSON overriding cohort/annotation settings")
    p.add_argument("--n-random", type=int, default=1000,
                   help="number of matched random datasets")
    p.add_argument("--seed", type=int, default=12345, help="random seed")
    p.add_argument("--no-umap", action="store_true",
                   help="ignore the Umap track (chromosome + length matching only)")
    return p.parse_args(argv)


def load_config(path):
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))  # deep copy
    if path:
        with open(path) as f:
            user = json.load(f)
        cfg.update(user)
    return cfg


def main(argv=None):
    args = parse_args(argv)
    cfg = load_config(args.config)

    out_dir = Path(args.out_dir)
    (out_dir / "source_data").mkdir(parents=True, exist_ok=True)
    ref_dir = Path(args.ref_dir)

    rng = random.Random(args.seed)
    np.random.seed(args.seed)

    chrom_sizes = load_chrom_sizes(ref_dir / cfg["chrom_sizes"])
    valid_chroms = set(chrom_sizes) & set(DEFAULT_HG38_CHROM_SIZES)

    annotations = {}
    for key, fname in cfg["annotation_files"].items():
        path = ref_dir / fname
        if not path.exists():
            log(f"[WARN] missing annotation BED: {path} -- '{key}' will be skipped")
            continue
        pairs = load_bed_pairs(path, valid_chroms)
        annotations[key] = IntervalSet.from_pairs(pairs)
        log(f"[LOAD] {key}: {sum(len(v) for v in pairs.values())} intervals")
    if not annotations:
        log("[ERROR] no annotation BED files could be loaded")
        sys.exit(1)
    if cfg["promoter_annotation"] not in annotations:
        log(f"[ERROR] promoter annotation '{cfg['promoter_annotation']}' not loaded; "
            f"set 'promoter_annotation' in --config to a loaded annotation")
        sys.exit(1)

    tss_arrays = load_protein_coding_tss(ref_dir / cfg["gencode_gtf"], valid_chroms)

    umap_pairs = None
    if not args.no_umap:
        umap_pairs = load_umap_pairs(ref_dir / cfg["umap_bedgraph"], valid_chroms)
    sampler = RandomIntervalSampler(chrom_sizes, umap_pairs, rng)

    numts = read_numt_calls(args.numt_root, cfg["cohorts"], valid_chroms)
    if numts.empty:
        log("[ERROR] no PASS_DISC_AND_SPLIT NUMT calls found under --numt-root")
        sys.exit(1)
    numts.to_csv(out_dir / "source_data" / "numts_used.tsv", sep="\t", index=False)
    log(f"[OK] total high-confidence NUMTs: {len(numts)} "
        f"across {numts['Cohort'].nunique()} cohorts")

    cohort_order = [c for c in cfg["cohort_order"] if c in set(numts["Cohort"])]

    enrich_df, tss_df, rand_overlap_long_df, rand_tss_long_df = run_analysis(
        numts, annotations, tss_arrays, sampler, args.n_random, cohort_order,
    )

    sd = out_dir / "source_data"
    key = cfg["promoter_annotation"]
    enrich_df.to_csv(sd / "annotation_enrichment.tsv", sep="\t", index=False)
    (enrich_df[enrich_df["Annotation"] == key]
     .to_csv(sd / "promoter_enrichment_by_cohort.tsv", sep="\t", index=False))
    (enrich_df[enrich_df["Annotation"] == key]
     [["Cohort", "N_NUMTs", "Observed_overlap", "Expected_overlap_mean",
       "Expected_overlap_sd", "Empirical_P", "Significance"]]
     .to_csv(sd / "promoter_observed_vs_expected.tsv", sep="\t", index=False))
    tss_df.to_csv(sd / "nearest_tss_distance.tsv", sep="\t", index=False)
    rand_overlap_long_df.to_csv(
        sd / "random_datasets_annotation_overlap_counts.tsv", sep="\t", index=False)
    rand_tss_long_df.to_csv(
        sd / "random_datasets_nearest_tss_per_dataset.tsv", sep="\t", index=False)

    make_figures(enrich_df, tss_df, cfg, out_dir)

    log("\n[DONE]")
    log(f"  figures:     {out_dir}")
    log(f"  source data: {sd}")
    sub = enrich_df[enrich_df["Annotation"] == key]
    log(f"\n{key} enrichment:")
    for _, r in sub.sort_values("Enrichment_fold", ascending=False).iterrows():
        log(f"  {r['Cohort']:<22} {r['Enrichment_fold']:.2f}-fold  "
            f"P={r['Empirical_P']:.3g} {r['Significance']}")


if __name__ == "__main__":
    main()
