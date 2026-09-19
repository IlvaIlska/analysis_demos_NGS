#!/usr/bin/env bash
#
# 01_align_wgbs_se.sh — QC, trim, and align a single-end WGBS run with Bismark, then deduplicate.
#
# The demo run (SRR1003257, GSE-linked SRP029721) is single-end bisulfite-seq. This produces a
# deduplicated BAM only — methylation calling is done downstream by dnmtools (stage 02), not by
# bismark_methylation_extractor. (The original 2023 script for this exact run used
# bismark_methylation_extractor with a --paired flag on single-end data, which was a bug; dnmtools
# format can consume the Bismark BAM directly, so that step is dropped here rather than copied.)
#
# INPUTS:
#   -1 <path>   Read 1 FASTQ (single-end)
#   -x <path>   Bismark genome folder (a directory prepared with `bismark_genome_preparation`,
#               e.g. .../mm10/ containing Bisulfite_Genome/)
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name; used as the prefix for every output file
#   -p <int>    Parallel instances for bismark (default: 4; bismark forks ~3-5 processes per
#               instance, so total core usage is roughly 3-5x this number)
#
# OUTPUTS (under -o):
#   trimmed/                              Trim Galore output (trimmed FASTQ + post-trim FastQC)
#   <name>_bismark_bt2.bam                Raw Bismark alignment
#   <name>_bismark_bt2.deduplicated.bam    Deduplicated alignment (final BAM used downstream)
#   *_bismark_bt2_SE_report.txt           Bismark alignment report
#   *.deduplication_report.txt            Bismark deduplication report
#
# REQUIRES: fastqc, trim_galore, bismark (with its bundled bowtie2), samtools

set -euo pipefail

parallel_instances=4

usage() { echo "Usage: $0 -1 <read1.fastq.gz> -x <bismark_genome_dir> -o <outdir> -s <sample_name> [-p <parallel_instances>]" >&2; exit 1; }

while getopts "1:x:o:s:p:" opt; do
  case "$opt" in
    1) read1=$OPTARG ;;
    x) genome_dir=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    p) parallel_instances=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${read1:-}" || -z "${genome_dir:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"
cd "$outdir"

# 1. Trim adapters + low-quality bases; Trim Galore also runs FastQC on the trimmed reads.
mkdir -p trimmed
trim_galore --gzip -o trimmed/ "$read1" --fastqc_args "-f fastq -o trimmed"

tRead1=$(ls trimmed/*trimmed.fq.gz)

# 2. Bismark alignment (single-end).
bismark --parallel "$parallel_instances" --nucleotide_coverage --prefix "$sample" "$genome_dir" --se "$tRead1"

bam=$(ls "${sample}"*bismark_bt2.bam)

# 3. Remove PCR duplicates.
deduplicate_bismark --bam "$bam"
