#!/usr/bin/env python3
"""
numtAnchorCluster.py

Clusters MT-nuclear discordant/split read pairs (produced by
NUMTs_detection_fixed.sh) into candidate NUMT insertion sites.

For each discordant read pair where one mate maps to the nuclear genome and
the other maps to MT, the nuclear-side coordinate is used as an "anchor".
Anchors on the same chromosome within --max-cluster-gap bp of each other are
grouped into a cluster. Each cluster is reported as a NUMT candidate if it is
supported by at least --min-disc-reads distinct read pairs, and is upgraded
to filter_status=PASS_DISC_AND_SPLIT if it also has >= --min-split-reads
split-read pairs corroborating the same region (otherwise
LOW_CONF_DISC_ONLY).

Usage:
  python3 numtAnchorCluster.py \\
      --sample SAMPLE --bam input.bam \\
      --disc sample.mt.disc.sam --split sample.mt.split.sam \\
      --out sample.NUMTs_candidates.tsv \\
      [--min-mapq 20] [--min-disc-reads 2] [--min-split-reads 1] \\
      [--max-cluster-gap 500]
"""

import argparse
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


MT_NAMES = {"chrM", "MT", "M"}


def is_primary_chrom(chrom: str) -> bool:
    if pd.isna(chrom):
        return False

    chrom = str(chrom)

    if chrom in MT_NAMES:
        return False

    if chrom in {"*", "="}:
        return False

    bad_patterns = [
        "random",
        "Un",
        "alt",
        "decoy",
        "HLA",
        "EBV",
        "GL",
        "KI",
        ".",
    ]

    for p in bad_patterns:
        if p in chrom:
            return False

    return True


def read_sam_as_df(path: str) -> pd.DataFrame:
    path = Path(path)

    rows = []
    with path.open() as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            rows.append(fields[:11])

    cols = [
        "QNAME", "FLAG", "RNAME", "POS", "MAPQ",
        "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"
    ]

    if len(rows) == 0:
        return pd.DataFrame(columns=cols)

    df = pd.DataFrame(rows, columns=cols)

    df["FLAG"] = pd.to_numeric(df["FLAG"], errors="coerce")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["PNEXT"] = pd.to_numeric(df["PNEXT"], errors="coerce")
    df["MAPQ"] = pd.to_numeric(df["MAPQ"], errors="coerce")
    df["READ_LEN"] = df["SEQ"].astype(str).str.len()

    return df


def cluster_positions(pos_list, max_gap: int):
    pos_list = sorted([int(x) for x in pos_list if pd.notna(x)])

    if len(pos_list) == 0:
        return []

    clusters = [[pos_list[0]]]

    for pos in pos_list[1:]:
        if pos - clusters[-1][-1] <= max_gap:
            clusters[-1].append(pos)
        else:
            clusters.append([pos])

    return clusters


def get_nuclear_anchor_rows(df: pd.DataFrame, min_mapq: int) -> pd.DataFrame:
    """
    Build nuclear-genome "anchor" coordinates from the discordant SAM.
    If the read itself maps to a nuclear chromosome, its POS is the anchor.
    If the read maps to MT and its mate maps to a nuclear chromosome, the
    mate's PNEXT is used as the anchor instead.
    """

    rows = []

    for _, r in df.iterrows():
        qname = r["QNAME"]
        rname = r["RNAME"]
        rnext = r["RNEXT"]
        pos = r["POS"]
        pnext = r["PNEXT"]
        mapq = r["MAPQ"]

        if pd.isna(mapq) or mapq < min_mapq:
            continue

        # case 1: read itself maps to nuclear genome, mate maps to MT
        if is_primary_chrom(rname) and rnext in MT_NAMES and pd.notna(pos):
            rows.append({
                "QNAME": qname,
                "chr": rname,
                "pos": int(pos),
                "MAPQ": float(mapq),
                "source": "RNAME_nuclear_RNEXT_MT"
            })

        # case 2: read itself maps to MT, mate maps to nuclear genome
        elif rname in MT_NAMES and is_primary_chrom(rnext) and pd.notna(pnext):
            rows.append({
                "QNAME": qname,
                "chr": rnext,
                "pos": int(pnext),
                "MAPQ": float(mapq),
                "source": "RNAME_MT_RNEXT_nuclear"
            })

    if len(rows) == 0:
        return pd.DataFrame(columns=["QNAME", "chr", "pos", "MAPQ", "source"])

    return pd.DataFrame(rows)


def get_mt_range(df_pair: pd.DataFrame):
    mt_df = df_pair[df_pair["RNAME"].isin(MT_NAMES)].copy()

    if mt_df.empty:
        return 0, 0

    mt_df = mt_df[pd.notna(mt_df["POS"])]

    if mt_df.empty:
        return 0, 0

    mt_start = int(mt_df["POS"].min())
    mt_end = int((mt_df["POS"] + mt_df["READ_LEN"] - 1).max())

    return mt_start, mt_end


def get_split_support_qnames(
    split_df: pd.DataFrame, chrom: str, start: int, end: int, min_mapq: int
) -> set:
    """
    Return the set of distinct QNAMEs in the split SAM that corroborate a
    candidate region. Returning the QNAME set (rather than just a count)
    lets the caller deduplicate against the discordant-read QNAMEs when
    computing total read support, since the same read pair can legitimately
    appear in both the discordant and split SAM files (e.g. one mate is
    part of a discordant pair *and* has a supplementary split alignment).
    """
    if split_df.empty:
        return set()

    sdf = split_df.copy()

    sdf = sdf[
        (sdf["MAPQ"].fillna(0) >= min_mapq) &
        (
            (
                (sdf["RNAME"] == chrom) &
                (sdf["POS"] >= start) &
                (sdf["POS"] <= end)
            )
            |
            (
                (sdf["RNEXT"] == chrom) &
                (sdf["PNEXT"] >= start) &
                (sdf["PNEXT"] <= end)
            )
        )
    ]

    if sdf.empty:
        return set()

    # Keep only split reads that actually involve MT (either mate on MT, or
    # a soft-/hard-clip in the CIGAR consistent with a split alignment).
    sdf = sdf[
        sdf["RNAME"].isin(MT_NAMES) |
        sdf["RNEXT"].isin(MT_NAMES) |
        sdf["CIGAR"].astype(str).str.contains("S|H", regex=True, na=False)
    ]

    return set(sdf["QNAME"].astype(str).unique())


def main():
    parser = argparse.ArgumentParser(
        description="Filter and cluster MT-nuclear discordant/split reads as NUMTs candidates."
    )

    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--disc", required=True)
    parser.add_argument("--split", required=True)
    parser.add_argument("--out", required=True)

    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-disc-reads", type=int, default=2)
    parser.add_argument("--min-split-reads", type=int, default=1)
    parser.add_argument("--max-cluster-gap", type=int, default=500)

    args = parser.parse_args()

    sample = args.sample
    bam = args.bam
    disc_path = Path(args.disc)
    split_path = Path(args.split)
    out_path = Path(args.out)

    logger.info("Sample: %s", sample)
    logger.info("Discordant SAM: %s", disc_path)
    logger.info("Split SAM: %s", split_path)
    logger.info("min MAPQ: %s", args.min_mapq)
    logger.info("min discordant reads: %s", args.min_disc_reads)
    logger.info("min split reads: %s", args.min_split_reads)

    disc_df = read_sam_as_df(disc_path)
    split_df = read_sam_as_df(split_path)

    out_cols = [
        "SampleID",
        "chr",
        "start",
        "end",
        "MT_start",
        "MT_end",
        "NUMT_discordant_reads",
        "NUMT_split_reads",
        "NUMT_total_reads",
        "mean_MAPQ",
        "filter_status",
        "discFile",
        "splitFile",
        "wgsBAM",
    ]

    if disc_df.empty:
        pd.DataFrame(columns=out_cols).to_csv(out_path, sep="\t", index=False)
        logger.info("No discordant reads. Empty output written.")
        return

    nuclear_df = get_nuclear_anchor_rows(disc_df, args.min_mapq)

    if nuclear_df.empty:
        pd.DataFrame(columns=out_cols).to_csv(out_path, sep="\t", index=False)
        logger.info("No nuclear anchors after MAPQ/chromosome filtering.")
        return

    results = []

    for chrom, group in nuclear_df.groupby("chr"):
        positions = group["pos"].dropna().astype(int).tolist()
        clusters = cluster_positions(positions, args.max_cluster_gap)

        for cluster_pos in clusters:
            raw_start = int(min(cluster_pos))
            raw_end = int(max(cluster_pos))

            cluster_df = group[
                (group["pos"] >= raw_start) &
                (group["pos"] <= raw_end)
            ]

            disc_qnames = set(cluster_df["QNAME"].astype(str))
            disc_support = len(disc_qnames)

            if disc_support < args.min_disc_reads:
                continue

            # Pad the cluster span so a downstream breakpoint search has
            # enough flanking sequence on both sides.
            start = max(0, raw_start - 500)
            end = raw_end + 650

            pair_df = disc_df[disc_df["QNAME"].astype(str).isin(disc_qnames)]
            mt_start, mt_end = get_mt_range(pair_df)

            split_qnames = get_split_support_qnames(
                split_df=split_df,
                chrom=chrom,
                start=start,
                end=end,
                min_mapq=args.min_mapq,
            )
            split_support = len(split_qnames)

            # A read pair can show up in both the discordant and split SAM
            # files (e.g. one mate is part of a discordant pair and also
            # has a supplementary split alignment). Dedup by QNAME so
            # NUMT_total_reads reflects distinct supporting read pairs
            # instead of double-counting them.
            total_support = len(disc_qnames | split_qnames)

            mean_mapq = round(float(cluster_df["MAPQ"].mean()), 2)

            if split_support >= args.min_split_reads:
                filter_status = "PASS_DISC_AND_SPLIT"
            else:
                filter_status = "LOW_CONF_DISC_ONLY"

            results.append([
                sample,
                chrom,
                start,
                end,
                mt_start,
                mt_end,
                disc_support,
                split_support,
                total_support,
                mean_mapq,
                filter_status,
                str(disc_path),
                str(split_path),
                str(bam),
            ])

    out_df = pd.DataFrame(results, columns=out_cols)

    if out_df.empty:
        out_df = pd.DataFrame(columns=out_cols)
    else:
        out_df = out_df.sort_values(
            ["filter_status", "chr", "start", "end"],
            ascending=[True, True, True, True]
        )

    out_df.to_csv(out_path, sep="\t", index=False)

    logger.info("DONE: %d candidates written to %s", len(out_df), out_path)
    if not out_df.empty:
        logger.info(
            "PASS_DISC_AND_SPLIT: %d",
            int((out_df["filter_status"] == "PASS_DISC_AND_SPLIT").sum())
        )
        logger.info(
            "LOW_CONF_DISC_ONLY: %d",
            int((out_df["filter_status"] == "LOW_CONF_DISC_ONLY").sum())
        )


if __name__ == "__main__":
    main()
