#!/usr/bin/env bash
#
# 04_derive_canyons.sh — split HMRs by length into "methylation canyons" (large hypomethylated
# domains) vs. shorter reference UMRs, generate length-matched random controls, and extract each
# region's boundary points for downstream profiling.
#
# RECONSTRUCTED THRESHOLDS: no script defining the canyon/UMR split survived — only the final
# .canyons and .ctrl.umr bed files did. The thresholds below were recovered empirically by
# measuring those surviving files' region-length distributions: canyons had min length 3502bp,
# ctrl.umr had max length 3499bp and min length 1000bp. That means the original split was almost
# certainly "HMR length >= 3500bp -> canyon; 1000bp <= length < 3500bp -> ctrl.umr; < 1000bp ->
# discarded" — which also matches the literature definition of a methylation canyon (Jeong et al.
# 2014, large hypomethylated domains >=3.5kb). These are documented, recovered constants, not
# arbitrary defaults.
#
# INPUTS:
#   -m <path>   HMR bed file (output of 03_call_hmrs.sh)
#   -g <path>   Genome chrom.sizes file (2-column: chrom, size) — needed for bedtools shuffle
#   -o <dir>    Output directory (created if absent)
#   -s <name>   Sample name; used as the prefix for every output file
#   -c <int>    Canyon minimum length in bp (default: 3500, recovered as above)
#   -u <int>    UMR minimum length in bp (default: 1000, recovered as above; HMRs shorter than
#               this are excluded from both sets)
#   -x <path>   Optional blacklist bed to exclude from the random shuffle (e.g. mm10 ENCODE
#               blacklist); omit to shuffle without exclusion
#
# OUTPUTS (under -o):
#   <name>.canyons / <name>.canyons.random           Canyons and their length-matched random control
#   <name>.ctrl.umr / <name>.ctrl.umr.random         Reference UMRs and their random control
#   start_<name>.canyons / end_<name>.canyons        1bp boundary points (and .random variants)
#   start_<name>.ctrl.umr / end_<name>.ctrl.umr      same, for the UMR set
#
# REQUIRES: bedtools

set -euo pipefail

canyon_min=3500
umr_min=1000
excl_args=()

usage() { echo "Usage: $0 -m <hmr.bed> -g <chrom.sizes> -o <outdir> -s <sample_name> [-c <canyon_min_bp>] [-u <umr_min_bp>] [-x <blacklist.bed>]" >&2; exit 1; }

while getopts "m:g:o:s:c:u:x:" opt; do
  case "$opt" in
    m) hmr=$OPTARG ;;
    g) chromsizes=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) sample=$OPTARG ;;
    c) canyon_min=$OPTARG ;;
    u) umr_min=$OPTARG ;;
    x) excl_args=(-excl "$OPTARG") ;;
    *) usage ;;
  esac
done

[[ -z "${hmr:-}" || -z "${chromsizes:-}" || -z "${outdir:-}" || -z "${sample:-}" ]] && usage

mkdir -p "$outdir"

canyons="${outdir}/${sample}.canyons"
umrs="${outdir}/${sample}.ctrl.umr"

# 1. Split by length.
awk -v min="$canyon_min" 'BEGIN{OFS="\t"} ($3-$2) >= min' "$hmr" > "$canyons"
awk -v lo="$umr_min" -v hi="$canyon_min" 'BEGIN{OFS="\t"} ($3-$2) >= lo && ($3-$2) < hi' "$hmr" > "$umrs"

# 2. Length-matched random controls (same interval lengths, shuffled genomic position).
bedtools shuffle -i "$canyons" -g "$chromsizes" -noOverlapping "${excl_args[@]}" > "${canyons}.random"
bedtools shuffle -i "$umrs" -g "$chromsizes" -noOverlapping "${excl_args[@]}" > "${umrs}.random"

# 3. Boundary extraction: 1bp points at each region's start and end, for both real and random sets.
for f in "$canyons" "${canyons}.random" "$umrs" "${umrs}.random"; do
  base=$(basename "$f")
  awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1}' "$f" > "${outdir}/start_${base}"
  awk 'BEGIN{OFS="\t"} {print $1, $3, $3+1}' "$f" > "${outdir}/end_${base}"
done
