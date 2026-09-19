#!/usr/bin/env bash
#
# 02_t4_spikein_normalization.sh — align the same trimmed reads to a T4 phage spike-in genome,
# use spike-in read depth to compute a per-sample scale factor, and apply it to the main-genome
# BAM from stage 01 to produce a spike-in-normalized bigWig (standard CUT&Tag/CUT&RUN
# normalization: spike-in DNA is added in constant amount per reaction, so its read depth
# estimates sample-to-sample efficiency differences that a simple RPKM normalization misses).
#
# INPUTS:
#   -1 <path>   Trimmed read 1 FASTQ (i.e. trimmed/*_R1*.gz from stage 01's output)
#   -2 <path>   Trimmed read 2 FASTQ
#   -x <path>   Bowtie2 index prefix for the spike-in genome (e.g. T4 phage; any T4 reference
#               works, such as RefSeq NC_000866.4)
#   -b <path>   Main-genome filtered+sorted BAM from stage 01 (<name>.s.filt.bam)
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name; used as the prefix for every output file
#   -p <int>    Threads (default: 4)
#
# OUTPUTS (under -o):
#   <name>_T4Map.rms.sam                  T4-aligned, sorted, deduplicated (has unmapped reads still)
#   <name>_T4Map.filt.bam                 Mapped-only spike-in reads (samtools view -F 0x04)
#   <name>_T4Map.bedpe / .filt.bedpe       Fragment pairs, filtered to same-chrom pairs < 1000bp
#   <name>_T4Map.frags.bedpe               Fragment start/end only, sorted
#   <name>.PicardMetrics.txt               Picard MarkDuplicates metrics (spike-in alignment)
#   <name>.mappedct / <name>.scale         Spike-in mapped read count and derived scale factor (10000 / count)
#   <name>.scaled.bw                       Main-genome bigWig, scaled by the spike-in factor
#
# REQUIRES: bowtie2, samtools, picard, bedtools, deeptools (bamCoverage)

set -euo pipefail

threads=4

usage() { echo "Usage: $0 -1 <trimmed_R1.gz> -2 <trimmed_R2.gz> -x <t4_bt2_index> -b <main_genome_filt.bam> -o <outdir> -s <sample_name> [-p <threads>]" >&2; exit 1; }

while getopts "1:2:x:b:o:s:p:" opt; do
  case "$opt" in
    1) read1=$OPTARG ;;
    2) read2=$OPTARG ;;
    x) t4_index=$OPTARG ;;
    b) main_bam=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    p) threads=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${read1:-}" || -z "${read2:-}" || -z "${t4_index:-}" || -z "${main_bam:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"
cd "$outdir"

# 1. Align trimmed reads to the T4 spike-in genome. Stringent end-to-end alignment with no
#    mixed/discordant pairs, since spike-in signal only needs concordant fragments.
bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant \
  --phred33 -I 10 -X 700 -p "$threads" -x "$t4_index" -1 "$read1" -2 "$read2" \
  -S "${sample}_T4Map.sam" &> "${sample}_T4Mapstats.txt"

# 2. Coordinate-sort and deduplicate.
picard SortSam I="${sample}_T4Map.sam" O="${sample}_T4Map.s.sam" SORT_ORDER=coordinate
picard MarkDuplicates \
  INPUT="${sample}_T4Map.s.sam" \
  OUTPUT="${sample}_T4Map.rms.sam" \
  METRICS_FILE="${sample}.PicardMetrics.txt" \
  REMOVE_DUPLICATES=true

# 3. Keep mapped reads only.
samtools view -bS -@ "$threads" -F 0x04 "${sample}_T4Map.rms.sam" > "${sample}_T4Map.filt.bam"

# 4. Fragment-level bedpe, filtered to same-chromosome pairs under 1kb.
bedtools bamtobed -i "${sample}_T4Map.filt.bam" -bedpe > "${sample}_T4Map.bedpe"
awk '$1==$4 && $6-$2 < 1000 {print $0}' "${sample}_T4Map.bedpe" > "${sample}_T4Map.filt.bedpe"
cut -f 1,2,6 "${sample}_T4Map.filt.bedpe" | sort -k1,1 -k2,2n -k3,3n > "${sample}_T4Map.frags.bedpe"

# 5. Scale factor from spike-in read depth: normalize every sample to 10000 spike-in reads.
samtools view -@ "$threads" -c "${sample}_T4Map.filt.bam" > "${sample}.mappedct"
awk '{print 10000/$1}' "${sample}.mappedct" > "${sample}.scale"
scale=$(cat "${sample}.scale")

# 6. Apply the scale factor to the main-genome BAM from stage 01.
samtools sort -@ "$threads" -o "${sample}_bt2.s.filt.bam" "$main_bam"
samtools index "${sample}_bt2.s.filt.bam"
bamCoverage -b "${sample}_bt2.s.filt.bam" -o "${sample}.scaled.bw" -p max --binSize 50 --scaleFactor "$scale"
