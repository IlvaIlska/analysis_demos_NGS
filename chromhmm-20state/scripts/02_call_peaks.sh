#!/usr/bin/env bash
#
# 02_call_peaks.sh — call narrowPeaks for one mark's single-end ChIP-seq BAM against the shared
# input control BAM. All 45 marks in this demo share one input (config/samples.tsv: SRR10192161).
#
# INPUTS:
#   -t <path>   Treatment (IP) BAM, final filtered+deduplicated BAM from stage 01
#   -c <path>   Control (input) BAM, same processing as -t (align the input run through stage 01 too)
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Mark name; used as the peak-name prefix and must match the mark name used in
#               config/cellmark_template.tab so stage 03 (BinarizeBed) can find it
#
# OUTPUTS (under -o):
#   <name>_nB_peaks.narrowPeak
#   <name>_nB_peaks.xls
#   <name>_nB_summits.bed
#   <name>_macs2Peak_summary.txt     MACS2 stderr log
#
# REQUIRES: macs2

set -euo pipefail

usage() { echo "Usage: $0 -t <treatment.bam> -c <control.bam> -o <outdir> -s <mark_name>" >&2; exit 1; }

while getopts "t:c:o:s:" opt; do
  case "$opt" in
    t) treatment=$OPTARG ;;
    c) control=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${treatment:-}" || -z "${control:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"

# -f BAM (not BAMPE): the demo data is single-end.
macs2 callpeak -t "$treatment" -c "$control" -g mm -f BAM \
  -n "${sample}_nB" --outdir "$outdir" -q 0.05 --keep-dup all \
  2> "${outdir}/${sample}_macs2Peak_summary.txt"
