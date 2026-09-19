#!/usr/bin/env bash
#
# 03_call_hmrs.sh — identify hypomethylated regions (HMRs) from a coverage-filtered, symmetric
# CpG methylation file, using dnmtools' 2-state HMM / beta-binomial emission model.
#
# INPUTS:
#   -m <path>   Coverage-filtered, symmetric .meth file (output of 02_build_methcounts.sh)
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name; used as the prefix for every output file
#
# OUTPUTS (under -o):
#   <name>.hmr           HMRs as a BED file: chrom, start, end, region label, #CpGs, strand (unused)
#   <name>.hmr.params    Learned HMM parameters (FG/BG alpha/beta, transition probabilities)
#   <name>.hmr.summary   Run summary (number of HMRs, mean size, etc.)
#
# REQUIRES: dnmtools

set -euo pipefail

usage() { echo "Usage: $0 -m <sym.10xCov.meth> -o <outdir> -s <sample_name>" >&2; exit 1; }

while getopts "m:o:s:" opt; do
  case "$opt" in
    m) methfile=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${methfile:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"

dnmtools hmr -v \
  -o "${outdir}/${sample}.hmr" \
  -p "${outdir}/${sample}.hmr.params" \
  -S "${outdir}/${sample}.hmr.summary" \
  "$methfile"
