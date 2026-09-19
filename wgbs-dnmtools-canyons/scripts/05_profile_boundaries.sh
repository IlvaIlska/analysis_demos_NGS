#!/usr/bin/env bash
#
# 05_profile_boundaries.sh — profile one or more signal bigWigs (5hmC, DNA methylation, etc.)
# across a window centered on canyon/UMR boundary points, and render a profile plot + heatmap.
#
# Typical use: point -s at the flagship pipeline's own sonicated-fraction 5hmC bigWig
# (../flagship-mnase-5hmc-diffbind/, stage 01/02 output) to see whether 5hmC is enriched or
# depleted right at methylation canyon edges — this is the cross-pipeline analysis this repo's
# root README alludes to for the canyon pipeline.
#
# A DNA methylation bigWig isn't produced by any earlier stage in this pipeline; build one from
# stage 02's .sym.meth file if you want to profile methylation itself alongside 5hmC, e.g.:
#   awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $5}' sample.sym.10xCov.meth | sort -k1,1 -k2,2n \
#     > sample.meth.bedGraph
#   bedGraphToBigWig sample.meth.bedGraph mm10.chrom.sizes sample.meth.bw
#
# INPUTS:
#   -r <path>       Boundary bed file (e.g. output of 04_derive_canyons.sh: start_<name>.canyons)
#   -s <path,...>   One or more signal bigWigs, comma-separated
#   -o <dir>        Output directory (created if absent)
#   -n <name>       Output file prefix
#   -d <int>        Flanking distance in bp, both upstream and downstream of the boundary
#                   (default: 3000)
#   -b <path>       Optional blacklist bed to exclude from the matrix (e.g. mm10 ENCODE blacklist)
#   -p <int>        Threads (default: 4)
#
# OUTPUTS (under -o):
#   <name>_dist<d>.mat            deepTools computeMatrix output
#   <name>_dist<d>_profile.pdf    Average signal profile across the window
#   <name>_dist<d>_heatmap.pdf    Per-region heatmap across the window
#
# REQUIRES: deeptools (computeMatrix, plotProfile, plotHeatmap)

set -euo pipefail

distance=3000
threads=4
blacklist_args=()

usage() { echo "Usage: $0 -r <boundary.bed> -s <signal1.bw,signal2.bw,...> -o <outdir> -n <name> [-d <flank_bp>] [-b <blacklist.bed>] [-p <threads>]" >&2; exit 1; }

while getopts "r:s:o:n:d:b:p:" opt; do
  case "$opt" in
    r) region=$OPTARG ;;
    s) signals=$OPTARG ;;
    o) outdir=$OPTARG ;;
    n) name=$OPTARG ;;
    d) distance=$OPTARG ;;
    b) blacklist_args=(-bl "$OPTARG") ;;
    p) threads=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${region:-}" || -z "${signals:-}" || -z "${outdir:-}" || -z "${name:-}" ]] && usage

mkdir -p "$outdir"
IFS=',' read -r -a signal_array <<< "$signals"

mat="${outdir}/${name}_dist${distance}.mat"

computeMatrix reference-point -R "$region" -S "${signal_array[@]}" \
  --referencePoint center -bs 10 --upstream "$distance" --downstream "$distance" \
  --sortRegions descend --smartLabels -p "$threads" --missingDataAsZero \
  "${blacklist_args[@]}" -o "$mat"

plotProfile -m "$mat" -o "${outdir}/${name}_dist${distance}_profile.pdf" --plotFileFormat pdf --yMin 0

plotHeatmap -m "$mat" -o "${outdir}/${name}_dist${distance}_heatmap.pdf" --plotFileFormat pdf \
  --colorMap viridis --whatToShow 'heatmap only' --zMin 0 --missingDataColor 0.5
