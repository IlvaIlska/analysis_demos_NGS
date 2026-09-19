#!/usr/bin/env bash
#
# 02_build_methcounts.sh — convert a deduplicated Bismark BAM into a coverage-filtered, symmetric
# CpG methylation-level file (dnmtools' standard methylome-construction chain: format -> sort ->
# uniq -> counts -> sym -> coverage filter).
#
# No original driver script for this conversion survived (only the final filtered .meth file did)
# — this follows dnmtools' own documented command chain (https://dnmtools.readthedocs.io/), not a
# recovered script.
#
# INPUTS:
#   -b <path>   Deduplicated Bismark BAM (output of 01_align_wgbs_se.sh)
#   -g <path>   Reference genome as a single FASTA file (dnmtools counts wants one file, not a
#               directory of per-chromosome FASTAs — if you only have the latter, build one with
#               `cat chr*.fa > mm10.fa`)
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name; used as the prefix for every output file
#   -c <int>    Minimum read coverage per CpG to keep (default: 10, matching the original run's
#               "10x coverage" filter documented in this pipeline's provenance notes)
#   -t <int>    Threads (default: 4)
#
# OUTPUTS (under -o):
#   <name>.formatted.bam         Bismark BAM reformatted to dnmtools' internal convention
#   <name>.formatted.sorted.bam  Coordinate-sorted version of the above
#   <name>.meth                  Per-cytosine (CpG-only) methylation levels (dnmtools counts)
#   <name>.sym.meth              Strand-symmetric CpG methylation levels (dnmtools sym)
#   <name>.sym.<c>xCov.meth      Final output: sym.meth filtered to >= -c coverage — this is the
#                                 file stage 03 (hmr calling) consumes
#
# REQUIRES: dnmtools (>=1.4, tested against 1.4.4), samtools

set -euo pipefail

mincov=10
threads=4

usage() { echo "Usage: $0 -b <deduplicated.bam> -g <genome.fa> -o <outdir> -s <sample_name> [-c <mincov>] [-t <threads>]" >&2; exit 1; }

while getopts "b:g:o:s:c:t:" opt; do
  case "$opt" in
    b) bam=$OPTARG ;;
    g) genome=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    c) mincov=$OPTARG ;;
    t) threads=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${bam:-}" || -z "${genome:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"
cd "$outdir"

# 1. Reformat the Bismark BAM into dnmtools' internal read convention (merges paired mates where
#    applicable; a no-op for single-end beyond format normalization).
dnmtools format -f bismark -single-end -t "$threads" -B "$bam" "${sample}.formatted.bam"

# 2. Coordinate-sort (dnmtools uniq requires reads sorted by chrom/start/end/strand).
samtools sort -@ "$threads" -o "${sample}.formatted.sorted.bam" "${sample}.formatted.bam"

# 3. Remove any remaining duplicate reads/fragments.
dnmtools uniq "${sample}.formatted.sorted.bam" "${sample}.uniq.bam"

# 4. Per-cytosine methylation levels, CpG context only.
dnmtools counts -t "$threads" -c "$genome" -n -o "${sample}.meth" "${sample}.uniq.bam"

# 5. Collapse to strand-symmetric CpG sites.
dnmtools sym -o "${sample}.sym.meth" "${sample}.meth"

# 6. Coverage filter.
awk -v c="$mincov" '$6 > c' "${sample}.sym.meth" > "${sample}.sym.${mincov}xCov.meth"
