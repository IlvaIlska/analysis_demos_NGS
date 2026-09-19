#!/usr/bin/env bash
#
# 01_trim_align_filter.sh — adapter-trim, align to the reference genome, restrict to properly
# paired reads, remove PCR/optical duplicates, and generate an RPKM-normalized bigWig.
#
# INPUTS:
#   -1 <path>   Read 1 FASTQ (gzipped)
#   -2 <path>   Read 2 FASTQ (gzipped)
#   -x <path>   Bowtie2 index prefix for the reference genome (e.g. .../mm10)
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name; used as the prefix for every output file
#   -p <int>    Threads (default: 8)
#
# OUTPUTS (under -o):
#   trimmed/                              Trim Galore output (trimmed FASTQs + post-trim FastQC)
#   <name>.s.bam(.bai)                    Sorted, unfiltered alignment
#   <name>.filt.bam                       Properly-paired reads only (samtools view -f 3)
#   <name>.s.filt.bam(.bai)               Sorted version of the above
#   <name>.s.rms.filt.bam(.bai)           Deduplicated (Picard MarkDuplicates), sorted, indexed — final BAM
#   <name>.PicardMetrics.txt              Picard MarkDuplicates metrics
#   <name>_INI_mapped_count.txt           samtools flagstat before dedup
#   <name>_FIN_mapped_count.txt           samtools flagstat after dedup
#   <name>.bw                             RPKM-normalized bigWig (deepTools bamCoverage)
#
# REQUIRES: trim_galore, bowtie2, samtools, picard (MarkDuplicates; $PICARD env var pointing at
#           the jar, or `picard` on PATH), deeptools (bamCoverage)

set -euo pipefail

threads=8

usage() { echo "Usage: $0 -1 <read1.fastq.gz> -2 <read2.fastq.gz> -x <bt2_index> -o <outdir> -s <sample_name> [-p <threads>]" >&2; exit 1; }

while getopts "1:2:x:o:s:p:" opt; do
  case "$opt" in
    1) read1=$OPTARG ;;
    2) read2=$OPTARG ;;
    x) bt2_index=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    p) threads=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${read1:-}" || -z "${read2:-}" || -z "${bt2_index:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"
cd "$outdir"

# 1. Trim adapters + low-quality bases; Trim Galore also runs FastQC on the trimmed reads.
mkdir -p trimmed
trim_galore --cores "$threads" --paired --gzip -o trimmed/ "$read1" "$read2" \
  --fastqc_args "-f fastq -o trimmed -t $threads"

tRead1=$(ls trimmed/*R1*.gz)
tRead2=$(ls trimmed/*R2*.gz)

# 2. Align to the reference genome.
bowtie2 --very-sensitive -p "$threads" -x "$bt2_index" -1 "$tRead1" -2 "$tRead2" \
  | samtools view -bS - > "${sample}.bam"
samtools sort -@ "$threads" -o "${sample}.s.bam" "${sample}.bam"
samtools index "${sample}.s.bam"

# 3. Restrict to properly-paired reads.
samtools view -bh -f 3 "${sample}.s.bam" > "${sample}.filt.bam"
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

# 6. RPKM-normalized bigWig for browser visualization.
bamCoverage -b "${sample}.s.rms.filt.bam" -o "${sample}.bw" --normalizeUsing RPKM -p max
