# WGBS/DNMTools 5hmC–5mC Methylation Canyon and HMR Calling

**Status: scaffolded, not yet built out.** The second of two secondary pipelines held in reserve
(see repo root README) — organized here so the directory structure and source material are
locked in, but the stage scripts haven't been cleaned up yet the way the
[flagship pipeline](../flagship-mnase-5hmc-diffbind/) has.

## What it does

Calls hypomethylated regions (HMRs) from WGBS methylation calls, then derives "methylation
canyons" — large hypomethylated domains bounded by sharp methylation transitions — and profiles
5hmC/nucleosome signal across canyon boundaries. Source:
`030325_callingMethCanyons_DNMTools_hmr/`.

## Pipeline stages (from source material)

1. **HMR calling**: `dnmtools hmr` on a per-CpG methylation-level file
   (`Casellas_nB_CpG_only.sym.10xCov.meth` — symmetric CpG calls, ≥10x coverage) → an HMR bed
   file (`Casellas_BS_naiveB.hmr`). Driver: `DNMTools_hmr.sh`.
2. **Canyon boundary derivation**: post-process the HMR calls with `bedtools` (see
   `bedtools_overlaps/`, `bedtools_opp_overlaps/`, `canyon_starts/`, `canyon_ends/`) to identify
   canyon boundaries — the source directory structure implies a "find HMR-to-HMR gaps above a
   size threshold, take the flanking transition points" approach, but the exact bedtools
   invocation needs to be re-derived from `beds/` before this is scripted cleanly.
3. **Boundary profiling**: `computeMatrix reference-point` (deepTools) centered on canyon
   boundaries, at several flank distances (1–5kb), against 5hmC/nucleosome bigWigs — produces
   profile plots and heatmaps (`canyon_scaleRegions.sh`).

## Demo data (candidate public accession — not yet wired up)

- [`SRP029721`](https://www.ncbi.nlm.nih.gov/sra/?term=SRP029721) — mouse B-cell resource
  dataset (~39 experiments: ChIA-PET, WGBS, H3K27Ac ChIP-seq, RNA-seq). This matches the
  "Casellas_nB_CpG_only" WGBS reference used throughout this pipeline — the WGBS runs within
  this accession are the direct public stand-in for `Casellas_nB_CpG_only.meth`.

The 5hmC/nucleosome bigWigs used for boundary profiling would come from the flagship pipeline's
own outputs (`../flagship-mnase-5hmc-diffbind/`), so this pipeline's demo is naturally downstream
of that one rather than fully independent.

## TODO before this is portfolio-ready

- [ ] Re-derive the exact canyon-boundary bedtools logic from `bedtools_overlaps/` /
      `bedtools_opp_overlaps/` output naming (not yet inspected in depth)
- [ ] Write a cleaned `01_call_hmrs.sh` (dnmtools hmr, parameterized on input .meth + output prefix)
- [ ] Write a cleaned `02_derive_canyons.sh` (the bedtools boundary logic above)
- [ ] Write a cleaned `03_profile_boundaries.sh` (deepTools computeMatrix/plotProfile/plotHeatmap,
      parameterized on flank distance instead of the original's fixed 1–5kb sweep)
- [ ] Confirm SRP029721's WGBS runs are processed the same way (symmetric CpG, ≥10x) as
      `Casellas_nB_CpG_only.sym.10xCov.meth` before treating them as a drop-in replacement
