#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_NUMTs_Fig_html.py

Aggregates *.NUMTs_candidates.tsv files produced by NUMTs_detection_fixed.sh
/ numtAnchorCluster.py (optionally across many samples, e.g. the
<output_root>/<sample>/ tree written by run_NUMTs_All_fixed.sh) and renders:

  - a circos-style static figure (PNG / SVG / PDF) linking mtDNA genes to
    the nuclear insertion sites, and
  - an interactive HTML viewer (Canvas-based) with per-sample and
    per-chromosome toggles.

Two aggregated views are produced by default: PASS_DISC_AND_SPLIT candidates
only, and PASS_DISC_AND_SPLIT + LOW_CONF_DISC_ONLY combined.

Usage:
  python3 run_NUMTs_Fig_html.py --input-dir <output_root> [--output-dir DIR]

--input-dir is searched recursively for files matching
"*.NUMTs_candidates.tsv" (the suffix used by NUMTs_detection_fixed.sh), so
it can point either at a single sample's output directory or at the root
directory containing one subdirectory per sample.
"""

import os
import glob
import math
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, PathPatch
from matplotlib.path import Path
from matplotlib import patheffects, colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams["font.family"] = ["DejaVu Sans", "Arial Unicode MS"]

# =========================================================
# CLI
# =========================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate NUMTs_candidates.tsv files produced by "
            "NUMTs_detection_fixed.sh / numtAnchorCluster.py and render "
            "circos-style figures (PNG/SVG/PDF) plus an interactive HTML "
            "viewer."
        )
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help=(
            "Root directory to search recursively for "
            "*.NUMTs_candidates.tsv files, e.g. the <output_root> passed "
            "to run_NUMTs_All_fixed.sh, or a single sample's output "
            "directory."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory under which the output_<label>/ figure/table "
            "folders are written (default: same as --input-dir)."
        ),
    )
    return parser.parse_args()


FILTER_SETS = {
    "PASS_DISC_AND_SPLIT": ["PASS_DISC_AND_SPLIT"],
    "PASS_PLUS_LOW_CONF_DISC_ONLY": [
        "PASS_DISC_AND_SPLIT",
        "LOW_CONF_DISC_ONLY",
    ],
}


# =========================================================
# Genome settings
# =========================================================

chr_len = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}

nuc_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

offsets = {}
cum_len = 0
for c in nuc_order:
    offsets[c] = cum_len
    cum_len += chr_len[c]

total_nuc_size = cum_len

ROTATE_DEG = 90.0
ROTATE_RAD = math.radians(ROTATE_DEG)


# =========================================================
# mtDNA annotation
# =========================================================

mt_genes = [
    ("Dloop", 1, 576),
    ("tRNA-Phe", 577, 647),
    ("RNR1", 648, 1601),
    ("tRNA-Val", 1602, 1670),
    ("RNR2", 1671, 3229),
    ("tRNA-Leu", 3230, 3304),
    ("ND1", 3307, 4262),
    ("tRNA-Ile", 4263, 4331),
    ("tRNA-Gln", 4329, 4400),
    ("tRNA-Met", 4402, 4469),
    ("ND2", 4470, 5511),
    ("tRNA-Trp", 5512, 5579),
    ("tRNA-Ala", 5587, 5655),
    ("tRNA-Asn", 5657, 5729),
    ("tRNA-Cys", 5761, 5826),
    ("tRNA-Tyr", 5826, 5891),
    ("CO1", 5904, 7445),
    ("tRNA-Ser", 7446, 7514),
    ("tRNA-Asp", 7518, 7585),
    ("COX2", 7586, 8269),
    ("tRNA-Lys", 8295, 8364),
    ("ATP8", 8366, 8572),
    ("ATP6", 8527, 9207),
    ("COX3", 9207, 9990),
    ("tRNA-Gly", 9991, 10058),
    ("ND3", 10059, 10404),
    ("tRNA-Arg", 10405, 10469),
    ("ND4L", 10470, 10766),
    ("ND4", 10760, 12137),
    ("tRNA-His", 12138, 12206),
    ("tRNA-Ser", 12207, 12265),
    ("tRNA-Leu", 12266, 12336),
    ("ND5", 12337, 14148),
    ("ND6", 14149, 14673),
    ("tRNA-Glu", 14674, 14742),
    ("CYTB", 14747, 15887),
    ("tRNA-Thr", 15888, 15953),
    ("tRNA-Pro", 15956, 16023),
    ("Dloop", 16024, 16569),
]


def mt_gene_from_pos(pos):
    pos = int(round(pos))
    for gene, s, e in mt_genes:
        if s <= pos <= e:
            return gene
    return "intergenic"


def pos_to_angle_mt(pos):
    pos = max(1, min(float(pos), chr_len["chrM"]))
    frac = (pos - 1) / (chr_len["chrM"] - 1)
    return math.radians(frac * 180.0) + ROTATE_RAD


def pos_to_angle_nuc(chr_name, pos):
    pos = max(0, min(float(pos), chr_len[chr_name]))
    frac = (offsets[chr_name] + pos) / total_nuc_size
    return math.radians(360.0 - frac * 180.0) + ROTATE_RAD


def normalize_chr(x):
    x = str(x).strip()

    if x in chr_len:
        return x

    if x in [str(i) for i in range(1, 23)]:
        return "chr" + x

    if x in ["X", "Y"]:
        return "chr" + x

    if x in ["M", "MT", "chrMT"]:
        return "chrM"

    return x


# =========================================================
# Colors
# =========================================================

def adjust_color(color, lighten=0.10, saturation_boost=1.08):
    rgb = np.array(mcolors.to_rgb(color))
    rgb = rgb + (1 - rgb) * lighten
    rgb = np.clip(rgb, 0, 1)
    hsv = mcolors.rgb_to_hsv(rgb)
    hsv[1] = min(1.0, hsv[1] * saturation_boost)
    return tuple(mcolors.hsv_to_rgb(hsv))


base_colors = plt.cm.Set3(np.linspace(0, 1, len(nuc_order)))
adj_colors = [adjust_color(c) for c in base_colors]
chr_colors = {c: adj_colors[i] for i, c in enumerate(nuc_order)}


def rgb_to_hex(rgb):
    return "#{:02x}{:02x}{:02x}".format(
        int(rgb[0] * 255),
        int(rgb[1] * 255),
        int(rgb[2] * 255),
    )


chr_colors_hex = {
    c: rgb_to_hex(mcolors.to_rgb(chr_colors[c]))
    for c in nuc_order
}


mt_gene_palette = {
    "Dloop": "#bdbdbd",
    "RNR1": "#80b1d3",
    "RNR2": "#80b1d3",
    "ND1": "#fb8072",
    "ND2": "#fb8072",
    "ND3": "#fb8072",
    "ND4": "#fb8072",
    "ND4L": "#fb8072",
    "ND5": "#fb8072",
    "ND6": "#fb8072",
    "CO1": "#8dd3c7",
    "COX2": "#8dd3c7",
    "COX3": "#8dd3c7",
    "CYTB": "#bebada",
    "ATP6": "#fdb462",
    "ATP8": "#fdb462",
}


def mt_gene_color(gene):
    if gene.startswith("tRNA"):
        return "#ffffb3"
    return mt_gene_palette.get(gene, "#d9d9d9")


# =========================================================
# Load function
# =========================================================

def load_data_by_filter(tsv_files, allowed_statuses):
    all_data = []

    required_cols = [
        "SampleID",
        "chr",
        "start",
        "end",
        "MT_start",
        "MT_end",
        "filter_status",
    ]

    for tsv_file in tsv_files:
        try:
            df = pd.read_csv(tsv_file, sep="\t")

            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                print(f"[SKIP] {os.path.basename(tsv_file)} missing {missing}")
                continue

            df = df[df["filter_status"].isin(allowed_statuses)].copy()

            if len(df) == 0:
                continue

            if "NUMT_reads" in df.columns:
                read_col = "NUMT_reads"

            elif "NUMT_total_reads" in df.columns:
                read_col = "NUMT_total_reads"

            elif (
                "NUMT_discordant_reads" in df.columns
                and "NUMT_split_reads" in df.columns
            ):
                df["NUMT_reads"] = (
                    pd.to_numeric(
                        df["NUMT_discordant_reads"],
                        errors="coerce",
                    ).fillna(0)
                    + pd.to_numeric(
                        df["NUMT_split_reads"],
                        errors="coerce",
                    ).fillna(0)
                )
                read_col = "NUMT_reads"

            else:
                print(f"[SKIP] {os.path.basename(tsv_file)} no NUMT read column")
                continue

            use_cols = [
                "SampleID",
                "chr",
                "start",
                "end",
                "MT_start",
                "MT_end",
                "filter_status",
                read_col,
            ]

            df = df[use_cols].copy()
            df = df.rename(columns={read_col: "NUMT_reads"})

            for c in ["start", "end", "MT_start", "MT_end", "NUMT_reads"]:
                df[c] = pd.to_numeric(df[c], errors="coerce")

            df = df.dropna(
                subset=[
                    "start",
                    "end",
                    "MT_start",
                    "MT_end",
                    "NUMT_reads",
                ]
            )

            df = df[df["NUMT_reads"] > 0].copy()

            if len(df) == 0:
                continue

            df["SampleID"] = df["SampleID"].astype(str)

            df.loc[df["SampleID"].isin(["", "nan", "None"]), "SampleID"] = (
                os.path.basename(tsv_file).replace(".tsv", "")
            )

            df["chr_mapped"] = df["chr"].apply(normalize_chr)
            df = df[df["chr_mapped"].isin(nuc_order)].copy()

            if len(df) == 0:
                continue

            all_data.append(df)

            print(
                f"[OK] {os.path.basename(tsv_file)} "
                f"rows={len(df)} "
                f"filters={allowed_statuses}"
            )

        except Exception as e:
            print(f"[ERROR] {os.path.basename(tsv_file)}")
            print(e)

    if len(all_data) == 0:
        return None

    return pd.concat(all_data, ignore_index=True)


# =========================================================
# Main plotting function
# =========================================================

def make_outputs(df_all, outdir, label_name):
    os.makedirs(outdir, exist_ok=True)

    grouped = (
        df_all.groupby(["chr_mapped", "MT_start"], as_index=False)
        .agg(
            NUMT_reads=("NUMT_reads", "sum"),
            start=("start", "mean"),
            end=("end", "mean"),
        )
    )

    grouped["mt_gene"] = grouped["MT_start"].apply(mt_gene_from_pos)

    max_reads = grouped["NUMT_reads"].max()
    n_samples = df_all["SampleID"].nunique()
    n_regions = len(grouped)
    total_reads = int(df_all["NUMT_reads"].sum())

    print(f"\n[MAKE] {label_name}")
    print(f"       samples: {n_samples}")
    print(f"       regions: {n_regions}")
    print(f"       total reads: {total_reads}")

    # =====================================================
    # TSV output
    # =====================================================

    rank_tsv = grouped.copy()
    rank_tsv["nuclear_midpoint"] = (
        (rank_tsv["start"] + rank_tsv["end"]) / 2
    ).round(1)

    rank_tsv["MT_start"] = rank_tsv["MT_start"].round().astype(int)

    rank_tsv = rank_tsv[
        [
            "mt_gene",
            "MT_start",
            "chr_mapped",
            "start",
            "end",
            "nuclear_midpoint",
            "NUMT_reads",
        ]
    ].rename(
        columns={
            "chr_mapped": "nuclear_chr",
            "start": "nuclear_start_mean",
            "end": "nuclear_end_mean",
            "NUMT_reads": "total_NUMT_reads",
        }
    )

    rank_tsv.sort_values(
        "total_NUMT_reads",
        ascending=False,
    ).to_csv(
        os.path.join(outdir, "NUMTs_top_links.tsv"),
        sep="\t",
        index=False,
    )

    rg = grouped.copy()
    rg["nuclear_midpoint"] = ((rg["start"] + rg["end"]) / 2).round(1)

    rg = (
        rg.groupby(["mt_gene", "chr_mapped"], as_index=False)
        .agg(
            total_NUMT_reads=("NUMT_reads", "sum"),
            link_count=("NUMT_reads", "size"),
            mean_nuclear_start=("start", "mean"),
            mean_nuclear_end=("end", "mean"),
            mean_nuclear_midpoint=("nuclear_midpoint", "mean"),
        )
        .sort_values(
            ["total_NUMT_reads", "link_count"],
            ascending=[False, False],
        )
        .reset_index(drop=True)
    )

    rg.rename(columns={"chr_mapped": "nuclear_chr"}).to_csv(
        os.path.join(outdir, "NUMTs_top_links_with_mt_gene.tsv"),
        sep="\t",
        index=False,
    )

    raw_out = df_all.copy()
    raw_out["mt_gene"] = raw_out["MT_start"].apply(mt_gene_from_pos)
    raw_out["nuclear_midpoint"] = ((raw_out["start"] + raw_out["end"]) / 2).round(1)

    raw_out.sort_values(
        "NUMT_reads",
        ascending=False,
    ).to_csv(
        os.path.join(outdir, "NUMTs_raw_records_sorted.tsv"),
        sep="\t",
        index=False,
    )

    # =====================================================
    # Static figure
    # =====================================================

    fig, ax = plt.subplots(
        figsize=(22, 22),
        subplot_kw=dict(aspect="equal"),
    )

    ax.set_facecolor("white")

    r_o = 1.12
    r_i = 0.94
    label_r = 1.26
    r_ctrl = 0.24

    for gene, s, e in mt_genes:
        a1 = math.degrees(pos_to_angle_mt(s))
        a2 = math.degrees(pos_to_angle_mt(e))

        ax.add_patch(
            Wedge(
                (0, 0),
                r_o,
                min(a1, a2),
                max(a1, a2),
                width=r_o - r_i,
                facecolor=mt_gene_color(gene),
                edgecolor="white",
                linewidth=1.8,
                alpha=0.98,
                zorder=3,
            )
        )

    for ch in nuc_order:
        s_ang = pos_to_angle_nuc(ch, 0)
        e_ang = pos_to_angle_nuc(ch, chr_len[ch])

        ax.add_patch(
            Wedge(
                (0, 0),
                r_o,
                min(math.degrees(s_ang), math.degrees(e_ang)),
                max(math.degrees(s_ang), math.degrees(e_ang)),
                width=r_o - r_i,
                facecolor="#d9d9d9",
                edgecolor="white",
                linewidth=2.3,
                alpha=0.95,
                zorder=3,
            )
        )

        mid = pos_to_angle_nuc(ch, chr_len[ch] / 2)

        ax.text(
            label_r * math.cos(mid),
            label_r * math.sin(mid),
            ch,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
            color="black",
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor="gray",
                alpha=0.95,
            ),
            zorder=20,
        )

    for gene, s, e in mt_genes:
        mid = (s + e) / 2
        ang = pos_to_angle_mt(mid)
        deg = math.degrees(ang) % 360

        rot = deg - 90
        ha = "left"

        if 90 < deg < 270:
            rot += 180
            ha = "right"

        ax.text(
            1.24 * math.cos(ang),
            1.24 * math.sin(ang),
            gene,
            ha=ha,
            va="center",
            fontsize=7.5 if gene.startswith("tRNA") else 9,
            fontweight="normal" if gene.startswith("tRNA") else "bold",
            rotation=rot,
            rotation_mode="anchor",
            zorder=15,
            bbox=dict(
                boxstyle="round,pad=0.08",
                facecolor="white",
                edgecolor="none",
                alpha=0.75,
            ),
        )

    mla = math.radians(90) + ROTATE_RAD

    ax.text(
        1.34 * math.cos(mla),
        1.34 * math.sin(mla),
        "chrM",
        ha="center",
        va="center",
        fontsize=22,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="gray",
            alpha=0.97,
        ),
        zorder=20,
    )

    for _, row in grouped.iterrows():
        mt_pos = row["MT_start"]
        nuc_pos = (row["start"] + row["end"]) / 2.0
        chr_name = row["chr_mapped"]

        mt_ang = pos_to_angle_mt(mt_pos)
        nuc_ang = pos_to_angle_nuc(chr_name, nuc_pos)

        x0 = r_i * math.cos(mt_ang)
        y0 = r_i * math.sin(mt_ang)

        x1 = r_ctrl * math.cos(mt_ang)
        y1 = r_ctrl * math.sin(mt_ang)

        x2 = r_ctrl * math.cos(nuc_ang)
        y2 = r_ctrl * math.sin(nuc_ang)

        x3 = r_i * math.cos(nuc_ang)
        y3 = r_i * math.sin(nuc_ang)

        read_ratio = row["NUMT_reads"] / max_reads if max_reads > 0 else 0

        lw = max(0.85, read_ratio * 14.5)
        alpha = min(0.95, 0.34 + read_ratio * 0.45)

        path = Path(
            [(x0, y0), (x1, y1), (x2, y2), (x3, y3)],
            [
                Path.MOVETO,
                Path.CURVE4,
                Path.CURVE4,
                Path.CURVE4,
            ],
        )

        ax.add_patch(
            PathPatch(
                path,
                facecolor="none",
                edgecolor=chr_colors[chr_name],
                lw=lw,
                alpha=alpha,
                capstyle="round",
                joinstyle="round",
                zorder=10,
                path_effects=[
                    patheffects.withStroke(
                        linewidth=lw + 0.7,
                        foreground="white",
                        alpha=0.30,
                    )
                ],
            )
        )

    plt.title(
        f"{label_name}\n"
        f"{n_regions} aggregated NUMT anchors from {n_samples} samples",
        fontsize=20,
        fontweight="bold",
        pad=35,
    )

    ax.set_xlim(-1.56, 1.56)
    ax.set_ylim(-1.56, 1.56)
    ax.axis("off")

    base = os.path.join(
        outdir,
        "NUMTs_semicircle_rotated90_segmented_chrM_grayRing_coloredLinks",
    )

    plt.savefig(base + ".png", dpi=450, bbox_inches="tight", facecolor="white")
    plt.savefig(base + ".svg", bbox_inches="tight", facecolor="white")

    with PdfPages(base + ".pdf") as pdf:
        pdf.savefig(plt.gcf(), dpi=350, bbox_inches="tight", facecolor="white")

    plt.close()

    # =====================================================
    # HTML output
    # =====================================================

    links_data = []

    for idx, row in grouped.iterrows():
        chr_name = row["chr_mapped"]
        mt_pos = row["MT_start"]
        nuc_pos = (row["start"] + row["end"]) / 2.0

        mt_angle = pos_to_angle_mt(mt_pos)
        nuc_angle = pos_to_angle_nuc(chr_name, nuc_pos)

        r_mt = 0.94
        r_nuc = 0.94
        r_ctrl = 0.24

        x0 = r_mt * math.cos(mt_angle)
        y0 = r_mt * math.sin(mt_angle)

        x1 = r_ctrl * math.cos(mt_angle)
        y1 = r_ctrl * math.sin(mt_angle)

        x2 = r_ctrl * math.cos(nuc_angle)
        y2 = r_ctrl * math.sin(nuc_angle)

        x3 = r_nuc * math.cos(nuc_angle)
        y3 = r_nuc * math.sin(nuc_angle)

        read_ratio = row["NUMT_reads"] / max_reads if max_reads > 0 else 0

        lw = max(0.85, read_ratio * 14.5)
        alpha = min(0.95, 0.34 + read_ratio * 0.45)

        sample_ids = (
            df_all[
                (df_all["chr_mapped"] == chr_name)
                & (df_all["MT_start"] == mt_pos)
            ]["SampleID"]
            .astype(str)
            .unique()
            .tolist()
        )

        links_data.append(
            {
                "id": f"link_{idx}",
                "p": [
                    round(x, 5)
                    for x in [
                        x0,
                        y0,
                        x1,
                        y1,
                        x2,
                        y2,
                        x3,
                        y3,
                    ]
                ],
                "c": chr_colors_hex[chr_name],
                "lw": round(lw, 2),
                "a": round(alpha, 3),
                "chr": chr_name,
                "g": row["mt_gene"],
                "mt": int(round(mt_pos)),
                "r": int(row["NUMT_reads"]),
                "s": sample_ids,
            }
        )

    mt_wedges = []

    for gene, s, e in mt_genes:
        mt_wedges.append(
            {
                "gene": gene,
                "t1": min(
                    math.degrees(pos_to_angle_mt(s)),
                    math.degrees(pos_to_angle_mt(e)),
                ),
                "t2": max(
                    math.degrees(pos_to_angle_mt(s)),
                    math.degrees(pos_to_angle_mt(e)),
                ),
                "c": mt_gene_color(gene),
            }
        )

    nuc_wedges = []

    for c in nuc_order:
        nuc_wedges.append(
            {
                "chr": c,
                "t1": min(
                    math.degrees(pos_to_angle_nuc(c, 0)),
                    math.degrees(pos_to_angle_nuc(c, chr_len[c])),
                ),
                "t2": max(
                    math.degrees(pos_to_angle_nuc(c, 0)),
                    math.degrees(pos_to_angle_nuc(c, chr_len[c])),
                ),
                "mid": math.degrees(pos_to_angle_nuc(c, chr_len[c] / 2)),
            }
        )

    samples_list = sorted(df_all["SampleID"].astype(str).unique().tolist())

    sample_totals = {
        str(k): int(v)
        for k, v in df_all.groupby("SampleID")["NUMT_reads"].sum().items()
    }

    chr_totals = {
        c: int(df_all[df_all["chr_mapped"] == c]["NUMT_reads"].sum())
        for c in nuc_order
    }

    config_json = json.dumps(
        {
            "links": links_data,
            "mt_wedges": mt_wedges,
            "nuc_wedges": nuc_wedges,
            "samples": samples_list,
            "chromosomes": nuc_order,
            "sample_totals": sample_totals,
            "chr_totals": chr_totals,
            "n_samples": n_samples,
            "n_regions": n_regions,
            "total_reads": total_reads,
            "chr_colors": chr_colors_hex,
            "label_name": label_name,
        },
        ensure_ascii=False,
    )

    html = f"""
<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="UTF-8">
<title>NUMTs Interactive Viewer - {label_name}</title>
<style>
body {{
    margin: 0;
    font-family: Arial, sans-serif;
    background: #111;
    color: #eee;
}}
header {{
    padding: 12px 20px;
    background: #1f1f1f;
    border-bottom: 1px solid #333;
}}
h1 {{
    margin: 0;
    font-size: 20px;
}}
#info {{
    margin-top: 6px;
    color: #aaa;
    font-size: 13px;
}}
#wrap {{
    display: flex;
    height: calc(100vh - 70px);
}}
#side {{
    width: 300px;
    overflow-y: auto;
    background: #181818;
    border-right: 1px solid #333;
    padding: 12px;
}}
#canvasWrap {{
    flex: 1;
    position: relative;
}}
canvas {{
    width: 100%;
    height: 100%;
    display: block;
}}
label {{
    display: block;
    font-size: 12px;
    margin: 4px 0;
}}
button {{
    margin: 3px;
    padding: 5px 8px;
}}
#tooltip {{
    position: fixed;
    display: none;
    background: rgba(0,0,0,0.85);
    color: white;
    padding: 8px 10px;
    border: 1px solid #777;
    border-radius: 6px;
    font-size: 12px;
    pointer-events: none;
    max-width: 260px;
}}
</style>
</head>
<body>

<header>
<h1>NUMTs Interactive Viewer - {label_name}</h1>
<div id="info">
Samples: {n_samples} | Aggregated regions: {n_regions} | Total reads: {total_reads}
</div>
</header>

<div id="wrap">
<div id="side">
<h3>Samples</h3>
<button onclick="selectAllSamples()">All</button>
<button onclick="clearSamples()">None</button>
<div id="sampleList"></div>

<hr>

<h3>Chromosomes</h3>
<button onclick="selectAllChroms()">All</button>
<button onclick="clearChroms()">None</button>
<div id="chromList"></div>
</div>

<div id="canvasWrap">
<canvas id="cv"></canvas>
</div>
</div>

<div id="tooltip"></div>

<script>
const C = {config_json};

const canvas = document.getElementById("cv");
const ctx = canvas.getContext("2d");
const tooltip = document.getElementById("tooltip");

let activeSamples = new Set(C.samples);
let activeChroms = new Set(C.chromosomes);
let hovered = null;

function resize() {{
    canvas.width = canvas.clientWidth * window.devicePixelRatio;
    canvas.height = canvas.clientHeight * window.devicePixelRatio;
    draw();
}}

window.addEventListener("resize", resize);

function center() {{
    return {{
        x: canvas.width / 2,
        y: canvas.height / 2
    }};
}}

function radius() {{
    return Math.min(canvas.width, canvas.height) * 0.38;
}}

function n2p(x, y) {{
    const c = center();
    const r = radius();
    return [
        c.x + x * r / 1.12,
        c.y - y * r / 1.12
    ];
}}

function hexToRgba(hex, alpha) {{
    const r = parseInt(hex.slice(1,3),16);
    const g = parseInt(hex.slice(3,5),16);
    const b = parseInt(hex.slice(5,7),16);
    return `rgba(${{r}},${{g}},${{b}},${{alpha}})`;
}}

function drawWedge(t1, t2, ro, ri, color) {{
    const c = center();
    const r = radius();
    const a1 = -t1 * Math.PI / 180;
    const a2 = -t2 * Math.PI / 180;

    ctx.beginPath();
    ctx.arc(c.x, c.y, r * ro / 1.12, a2, a1, false);
    ctx.arc(c.x, c.y, r * ri / 1.12, a1, a2, true);
    ctx.closePath();
    ctx.fillStyle = color;
    ctx.fill();
    ctx.strokeStyle = "#fff";
    ctx.lineWidth = 1;
    ctx.stroke();
}}

function draw() {{
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = "#111";
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    for (const w of C.mt_wedges) {{
        drawWedge(w.t1, w.t2, 1.12, 0.94, w.c);
    }}

    for (const w of C.nuc_wedges) {{
        drawWedge(w.t1, w.t2, 1.12, 0.94, "#d9d9d9");
    }}

    ctx.lineCap = "round";
    ctx.lineJoin = "round";

    for (const link of C.links) {{
        const sampleOK = link.s.some(s => activeSamples.has(s));
        const chrOK = activeChroms.has(link.chr);

        if (!sampleOK || !chrOK) continue;

        const p = link.p;
        const p0 = n2p(p[0], p[1]);
        const p1 = n2p(p[2], p[3]);
        const p2 = n2p(p[4], p[5]);
        const p3 = n2p(p[6], p[7]);

        ctx.beginPath();
        ctx.moveTo(p0[0], p0[1]);
        ctx.bezierCurveTo(
            p1[0], p1[1],
            p2[0], p2[1],
            p3[0], p3[1]
        );

        ctx.strokeStyle = link.id === hovered
            ? "#ffffff"
            : hexToRgba(link.c, link.a);

        ctx.lineWidth = link.id === hovered
            ? link.lw * 2.2
            : link.lw;

        ctx.stroke();
    }}

    drawLabels();
}}

function drawLabels() {{
    const c = center();
    const r = radius();

    ctx.fillStyle = "#eee";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.font = "bold 12px Arial";

    for (const w of C.nuc_wedges) {{
        const a = -w.mid * Math.PI / 180;
        const x = c.x + Math.cos(a) * r * 1.26 / 1.12;
        const y = c.y + Math.sin(a) * r * 1.26 / 1.12;
        ctx.fillText(w.chr, x, y);
    }}

    ctx.font = "bold 18px Arial";
    ctx.fillText("chrM", c.x, c.y - r * 1.20);
}}

function distPointToSegment(px, py, x1, y1, x2, y2) {{
    const A = px - x1;
    const B = py - y1;
    const Cc = x2 - x1;
    const D = y2 - y1;

    const dot = A * Cc + B * D;
    const lenSq = Cc * Cc + D * D;

    let param = -1;

    if (lenSq !== 0) param = dot / lenSq;

    let xx, yy;

    if (param < 0) {{
        xx = x1;
        yy = y1;
    }} else if (param > 1) {{
        xx = x2;
        yy = y2;
    }} else {{
        xx = x1 + param * Cc;
        yy = y1 + param * D;
    }}

    const dx = px - xx;
    const dy = py - yy;

    return Math.sqrt(dx * dx + dy * dy);
}}

canvas.addEventListener("mousemove", e => {{
    const rect = canvas.getBoundingClientRect();

    const mx = (e.clientX - rect.left) * window.devicePixelRatio;
    const my = (e.clientY - rect.top) * window.devicePixelRatio;

    let hit = null;
    let best = 999999;

    for (const link of C.links) {{
        const sampleOK = link.s.some(s => activeSamples.has(s));
        const chrOK = activeChroms.has(link.chr);

        if (!sampleOK || !chrOK) continue;

        const p = link.p;
        const p0 = n2p(p[0], p[1]);
        const p3 = n2p(p[6], p[7]);

        const d = distPointToSegment(
            mx,
            my,
            p0[0],
            p0[1],
            p3[0],
            p3[1]
        );

        if (d < best && d < 12) {{
            best = d;
            hit = link;
        }}
    }}

    if (hit) {{
        hovered = hit.id;

        tooltip.style.display = "block";
        tooltip.style.left = e.clientX + 15 + "px";
        tooltip.style.top = e.clientY + 15 + "px";

        tooltip.innerHTML = `
            <b>${{hit.g}}</b><br>
            mt position: ${{hit.mt}}<br>
            nuclear chr: ${{hit.chr}}<br>
            reads: ${{hit.r}}<br>
            samples: ${{hit.s.join(", ")}}
        `;
    }} else {{
        hovered = null;
        tooltip.style.display = "none";
    }}

    draw();
}});

canvas.addEventListener("mouseleave", () => {{
    hovered = null;
    tooltip.style.display = "none";
    draw();
}});

function buildSidebars() {{
    const sampleList = document.getElementById("sampleList");
    sampleList.innerHTML = "";

    const sampleEntries = Object.entries(C.sample_totals)
        .sort((a, b) => b[1] - a[1]);

    for (const [sample, reads] of sampleEntries) {{
        const label = document.createElement("label");
        label.innerHTML = `
            <input type="checkbox" checked data-sample="${{sample}}">
            ${{sample}} (${{reads}})
        `;

        label.querySelector("input").addEventListener("change", e => {{
            if (e.target.checked) activeSamples.add(sample);
            else activeSamples.delete(sample);
            draw();
        }});

        sampleList.appendChild(label);
    }}

    const chromList = document.getElementById("chromList");
    chromList.innerHTML = "";

    for (const chr of C.chromosomes) {{
        const reads = C.chr_totals[chr] || 0;
        const color = C.chr_colors[chr];

        const label = document.createElement("label");
        label.innerHTML = `
            <input type="checkbox" checked data-chr="${{chr}}">
            <span style="color:${{color}}">■</span>
            ${{chr}} (${{reads}})
        `;

        label.querySelector("input").addEventListener("change", e => {{
            if (e.target.checked) activeChroms.add(chr);
            else activeChroms.delete(chr);
            draw();
        }});

        chromList.appendChild(label);
    }}
}}

function selectAllSamples() {{
    activeSamples = new Set(C.samples);
    document.querySelectorAll("[data-sample]").forEach(x => x.checked = true);
    draw();
}}

function clearSamples() {{
    activeSamples = new Set();
    document.querySelectorAll("[data-sample]").forEach(x => x.checked = false);
    draw();
}}

function selectAllChroms() {{
    activeChroms = new Set(C.chromosomes);
    document.querySelectorAll("[data-chr]").forEach(x => x.checked = true);
    draw();
}}

function clearChroms() {{
    activeChroms = new Set();
    document.querySelectorAll("[data-chr]").forEach(x => x.checked = false);
    draw();
}}

buildSidebars();
resize();
</script>

</body>
</html>
"""

    with open(
        os.path.join(outdir, "NUMTs_interactive_viewer.html"),
        "w",
        encoding="utf-8",
    ) as f:
        f.write(html)

    print(f"[DONE] {label_name}")
    print(f"       output: {outdir}")


# =========================================================
# Main
# =========================================================

def main():
    args = parse_args()

    base_dir = args.input_dir
    out_root = args.output_dir if args.output_dir else args.input_dir

    tsv_files = sorted(
        glob.glob(
            os.path.join(base_dir, "**", "*.NUMTs_candidates.tsv"),
            recursive=True,
        )
    )

    print(f"[INFO] input dir: {base_dir}")
    print(f"[INFO] TSV files found: {len(tsv_files)}")

    if len(tsv_files) == 0:
        raise FileNotFoundError(
            f"No *.NUMTs_candidates.tsv files found under {base_dir}"
        )

    os.makedirs(out_root, exist_ok=True)

    for label_name, statuses in FILTER_SETS.items():
        print("\n====================================================")
        print(f"[RUN] {label_name}")
        print(f"[FILTER] {statuses}")
        print("====================================================")

        df = load_data_by_filter(tsv_files, statuses)

        if df is None or len(df) == 0:
            print(f"[SKIP] No valid rows for {label_name}")
            continue

        outdir = os.path.join(out_root, f"output_{label_name}")

        make_outputs(
            df_all=df,
            outdir=outdir,
            label_name=label_name,
        )

    print("\nALL DONE")


if __name__ == "__main__":
    main()
