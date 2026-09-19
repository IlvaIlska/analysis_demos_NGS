#!/usr/bin/env bash
#
# 01_align_single_end.sh — adapter-trim, align a single-end ChIP-seq run to the reference
# genome, filter to mapped reads, and remove PCR/optical duplicates.
#
# The demo data (GSE82144/SRP075985) is single-end HiSeq 2500 ChIP-seq, unlike the flagship
# pipeline's paired-end CUT&Tag — hence a separate alignment script rather than reusing
# ../flagship-mnase-5hmc-diffbind/scripts/01_trim_align_filter.sh.
#
# INPUTS:
#   -1 <path>   Read 1 FASTQ (gzipped; single-end, so only one file)
#   -x <path>   Bowtie2 index prefix for the reference genome (e.g. .../mm10)
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name (e.g. mark name); used as the prefix for every output file
#   -p <int>    Threads (default: 8)
#
# OUTPUTS (under -o):
#   trimmed/                        Trim Galore output (trimmed FASTQ + post-trim FastQC)
#   <name>.s.bam(.bai)              Sorted, unfiltered alignment
#   <name>.rms.filt.bam             Mapped reads only, deduplicated (final BAM used for peak calling)
#   <name>.s.rms.filt.bam(.bai)     Sorted version of the above
#   <name>.PicardMetrics.txt        Picard MarkDuplicates metrics
#   <name>_INI_mapped_count.txt     samtools flagstat before filtering
#   <name>_FIN_mapped_count.txt     samtools flagstat after filtering+dedup
#
# REQUIRES: trim_galore, bowtie2, samtools, picard (MarkDuplicates)

set -euo pipefail

threads=8

usage() { echo "Usage: $0 -1 <read1.fastq.gz> -x <bt2_index> -o <outdir> -s <sample_name> [-p <threads>]" >&2; exit 1; }

while getopts "1:x:o:s:p:" opt; do
  case "$opt" in
    1) read1=$OPTARG ;;
    x) bt2_index=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    p) threads=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${read1:-}" || -z "${bt2_index:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"
cd "$outdir"

# 1. Trim adapters + low-quality bases; Trim Galore also runs FastQC on the trimmed reads.
mkdir -p trimmed
trim_galore --cores "$threads" --gzip -o trimmed/ "$read1" \
  --fastqc_args "-f fastq -o trimmed -t $threads"

tRead1=$(ls trimmed/*trimmed.fq.gz)

# 2. Align to the reference genome (single-end).
bowtie2 --very-sensitive -p "$threads" -x "$bt2_index" -U "$tRead1" \
  | samtools view -bS - > "${sample}.bam"
samtools sort -@ "$threads" -o "${sample}.s.bam" "${sample}.bam"
samtools index "${sample}.s.bam"

# 3. Keep mapped reads only (-F 4).
samtools view -bh -F 4 "${sample}.s.bam" > "${sample}.filt.bam"
samtools sort -@ "$threads" -o "${sample}.s.filt.bam" "${sample}.filt.bam"
samtools index "${sample}.s.filt.bam"

# 4. Remove PCR/optical duplicates.
picard MarkDuplicates \
  INPUT="${sample}.s.filt.bam" \
  OUTPUT="${sample}.rms.filt.bam" \
  METRICS_FILE="${sample}.PicardMetrics.txt" \
  REMOVE_DUPLICATES=true
samtools sort -@ "$threads" -o "${sample}.s.rms.filt.bam" "${sample}.rms.filt.bam"
samtools index "${sample}.s.rms.filt.bam"

# 5. Mapping stats before/after filtering, for QC.
samtools flagstat "${sample}.s.bam" > "${sample}_INI_mapped_count.txt"
samtools flagstat "${sample}.s.rms.filt.bam" > "${sample}_FIN_mapped_count.txt"
