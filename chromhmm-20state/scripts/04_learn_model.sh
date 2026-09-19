#!/usr/bin/env bash
#
# 04_learn_model.sh — fit a ChromHMM chromatin-state model from the binarized mark matrix and
# produce the full state-calling output (segmentation, emission/transition parameters, TSS/TES
# enrichment plots, browser/dense bed tracks, summary webpage).
#
# INPUTS:
#   -i <dir>    Binarized input directory (output of 03_binarize.sh)
#   -o <dir>    Output directory (created if absent)
#   -n <int>    Number of states (default: 20, matching the original run)
#   -a <str>    Genome assembly name known to ChromHMM (default: mm10 — ChromHMM ships mm10
#               chrom sizes/coordinates/anchor files, so no extra reference download is needed)
#   -p <int>    Max parallel processors (default: 4)
#   -j <path>   Path to ChromHMM.jar (default: $CHROMHMM_JAR)
#
# OUTPUTS (under -o):
#   model_<n>.txt                     Learned HMM model (emission/transition parameters)
#   emissions_<n>.txt/.png/.svg       Emission parameter table + heatmap
#   transitions_<n>.txt/.png/.svg     Transition parameter table + heatmap
#   nB_<n>_segments.bed.gz            Per-bin state calls (genome-wide segmentation)
#   nB_<n>_dense.bed.gz               UCSC-browser-ready dense state track
#   nB_<n>_RefSeqTSS/TES_neighborhood.txt/.png/.svg   State enrichment around TSS/TES
#   nB_<n>_overlap.txt/.png/.svg      State overlap enrichment against ChromHMM's built-in annotations
#   webpage_<n>.html                  Summary report tying all of the above together
#   POSTERIOR/                        Per-chromosome posterior state probabilities
#
# REQUIRES: java, ChromHMM.jar (LearnModel mode)

set -euo pipefail

numstates=20
assembly=mm10
maxprocessors=4
jar="${CHROMHMM_JAR:-}"

usage() { echo "Usage: $0 -i <binarized_dir> -o <outdir> [-n <numstates>] [-a <assembly>] [-p <maxprocessors>] [-j <ChromHMM.jar>]" >&2; exit 1; }

while getopts "i:o:n:a:p:j:" opt; do
  case "$opt" in
    i) inputdir=$OPTARG ;;
    o) outdir=$OPTARG ;;
    n) numstates=$OPTARG ;;
    a) assembly=$OPTARG ;;
    p) maxprocessors=$OPTARG ;;
    j) jar=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${inputdir:-}" || -z "${outdir:-}" ]] && usage
[[ -z "$jar" ]] && { echo "No ChromHMM.jar path given (-j) and \$CHROMHMM_JAR is unset." >&2; usage; }

mkdir -p "$outdir"

java -jar "$jar" LearnModel -p "$maxprocessors" \
  "$inputdir" "$outdir" "$numstates" "$assembly"
