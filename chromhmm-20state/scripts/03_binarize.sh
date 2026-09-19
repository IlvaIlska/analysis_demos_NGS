#!/usr/bin/env bash
#
# 03_binarize.sh — convert per-mark narrowPeak calls into ChromHMM's binary presence/absence
# matrix (BinarizeBed), keyed by a cellmark table.
#
# INPUTS:
#   -i <dir>    Input directory containing all narrowPeak files named exactly as listed in the
#               cellmark table's 3rd column (e.g. config/cellmark_template.tab)
#   -c <path>   Cellmark table: <celltype> <mark> <narrowPeak filename>, one row per mark
#               (config/cellmark_template.tab)
#   -o <dir>    Output directory for binarized files (created if absent)
#   -g <path>   Chromosome length file (2-column: chrom, size) — this repo doesn't ship one;
#               use ChromHMM's bundled one, e.g.
#               /fs/ess/PAS1713/0-common/1-softwares/ChromHMM/CHROMSIZES/mm10.txt, or point at
#               your own mm10.chrom.sizes
#   -b <int>    Bin size in bp (default: 200, ChromHMM's own default)
#   -j <path>   Path to ChromHMM.jar (default: assumes `chromhmm.jar` on PATH via $CHROMHMM_JAR)
#
# OUTPUTS (under -o):
#   nB_<chrom>_binary.txt.gz   One file per chromosome: header (celltype, chrom; mark names),
#                              then one row per bin, 0/1 per mark
#
# REQUIRES: java, ChromHMM.jar (BinarizeBed mode)

set -euo pipefail

binsize=200
jar="${CHROMHMM_JAR:-}"

usage() { echo "Usage: $0 -i <narrowpeak_dir> -c <cellmark_table> -o <outdir> -g <chrom_sizes> [-b <binsize>] [-j <ChromHMM.jar>]" >&2; exit 1; }

while getopts "i:c:o:g:b:j:" opt; do
  case "$opt" in
    i) inputbeddir=$OPTARG ;;
    c) cellmarktable=$OPTARG ;;
    o) outdir=$OPTARG ;;
    g) chromsizes=$OPTARG ;;
    b) binsize=$OPTARG ;;
    j) jar=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${inputbeddir:-}" || -z "${cellmarktable:-}" || -z "${outdir:-}" || -z "${chromsizes:-}" ]] && usage
[[ -z "$jar" ]] && { echo "No ChromHMM.jar path given (-j) and \$CHROMHMM_JAR is unset." >&2; usage; }

mkdir -p "$outdir"

java -jar "$jar" BinarizeBed -b "$binsize" -peaks \
  "$chromsizes" "$inputbeddir" "$cellmarktable" "$outdir"
