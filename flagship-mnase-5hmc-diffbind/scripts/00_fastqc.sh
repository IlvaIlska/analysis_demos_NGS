#!/usr/bin/env bash
#
# 00_fastqc.sh — raw-read quality control.
#
# INPUTS:
#   -1 <path>   Read 1 FASTQ (gzipped)
#   -2 <path>   Read 2 FASTQ (gzipped)
#   -o <dir>    Output directory (created if absent)
#   -p <int>    Threads (default: 4)
#
# OUTPUTS (under -o):
#   <read1_basename>_fastqc.html / .zip
#   <read2_basename>_fastqc.html / .zip
#
# REQUIRES: fastqc
#
# Example:
#   ./00_fastqc.sh -1 submono_1590_R1.fastq.gz -2 submono_1590_R2.fastq.gz \
#       -o results/submono_1590/qc_raw

set -euo pipefail

threads=4

usage() { echo "Usage: $0 -1 <read1.fastq.gz> -2 <read2.fastq.gz> -o <outdir> [-p <threads>]" >&2; exit 1; }

while getopts "1:2:o:p:" opt; do
  case "$opt" in
    1) read1=$OPTARG ;;
    2) read2=$OPTARG ;;
    o) outdir=$OPTARG ;;
    p) threads=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${read1:-}" || -z "${read2:-}" || -z "${outdir:-}" ]] && usage

mkdir -p "$outdir"
fastqc -f fastq -o "$outdir" --threads "$threads" "$read1" "$read2"
