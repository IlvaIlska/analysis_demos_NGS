#!/usr/bin/env bash
#
# 03_macs2_call_peaks.sh — call 5hmC peaks for one IP sample against its matched input control.
#
# INPUTS:
#   -t <path>   Treatment (IP) BAM, final filtered+deduplicated BAM from stage 01
#   -c <path>   Control (input) BAM, same processing as -t
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name; used as the peak-name prefix
#
# OUTPUTS (under -o):
#   <name>_macs2_q0.05_peaks.narrowPeak
#   <name>_macs2_q0.05_peaks.xls
#   <name>_macs2_q0.05_summits.bed
#   <name>_macs2Peak_summary.txt          MACS2 stderr log
#
# REQUIRES: macs2

set -euo pipefail

usage() { echo "Usage: $0 -t <treatment.bam> -c <control.bam> -o <outdir> -s <sample_name>" >&2; exit 1; }

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

macs2 callpeak -t "$treatment" -c "$control" -g mm -f BAMPE \
  -n "${sample}_macs2_q0.05" --outdir "$outdir" -q 0.05 --keep-dup all \
  2> "${outdir}/${sample}_macs2Peak_summary.txt"
