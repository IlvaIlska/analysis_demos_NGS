# WGBS/DNMTools 5hmC–5mC Methylation Canyon and HMR Calling

Calls hypomethylated regions (HMRs) from WGBS methylation calls in naïve mouse B cells, splits
them into large "methylation canyons" vs. shorter reference UMRs, generates length-matched random
controls, and profiles 5hmC/methylation signal across canyon boundaries. Source:
`030325_callingMethCanyons_DNMTools_hmr/` (see repo root README for provenance).

## Pipeline stages

| Stage | Script | What it does |
|---|---|---|
| 1 | `scripts/01_align_wgbs_se.sh` | Trim (Trim Galore) → align single-end WGBS to mm10 (Bismark) → dedup |
| 2 | `scripts/02_build_methcounts.sh` | dnmtools `format` → sort → `uniq` → `counts` → `sym` → coverage filter (≥10x) |
| 3 | `scripts/03_call_hmrs.sh` | dnmtools `hmr`: 2-state HMM hypomethylated-region calling |
| 4 | `scripts/04_derive_canyons.sh` | Split HMRs by length into canyons (≥3.5kb) vs. reference UMRs (1–3.5kb) → length-matched random controls (`bedtools shuffle`) → 1bp boundary extraction |
| 5 | `scripts/05_profile_boundaries.sh` | deepTools `computeMatrix reference-point`/`plotProfile`/`plotHeatmap` centered on canyon/UMR boundaries |

Each script has a header block listing its declared inputs, outputs, and parameters.

## Two reconstructions, one recovered from evidence in the surviving data

**Stage 2** (the Bismark BAM → dnmtools methylation-counts chain) has no surviving driver script
— only the final filtered `.meth` file did. `02_build_methcounts.sh` follows dnmtools' own
documented command chain (`format` → sort → `uniq` → `counts` → `sym`), not a recovered script.
The original 2023-era script for this exact run (naive-B WGBS, SRR1003257) used
`bismark_methylation_extractor --paired` on single-end data, which — paired flag on single-end
input — looks like a bug; this pipeline uses the modern dnmtools chain instead of reproducing it.

**Stage 4** (the canyon/UMR length split) also has no surviving driver script — but *does* have
surviving output (`Casellas_BS_naiveB.canyons` and `Casellas_BS_naiveB.ctrl.umr`), which let the
threshold be recovered empirically rather than guessed: measuring those files' region lengths
gives canyons min length **3502bp** and ctrl.umr max length **3499bp** / min length **1000bp** —
i.e. the original split was almost certainly "≥3500bp → canyon, 1000–3499bp → reference UMR,
<1000bp → discarded." This also matches the literature definition of a methylation canyon (Jeong
et al. 2014, large hypomethylated domains ≥3.5kb), so it's used here as a documented, recovered
constant rather than an arbitrary default. Verified against real `bedtools` with a synthetic
4-region HMR file before committing (correct 3-way split, length-preserving shuffle, correct
1bp boundary points).

## Demo data

**Public accession:** [`SRP029721`](https://www.ncbi.nlm.nih.gov/sra/?term=SRP029721) (Casellas
lab mouse B-cell WGBS resource) — specifically **`SRR1003257`**, "wild-type resting B cells,"
single-end Illumina HiSeq 2000. This is the exact run the original `Casellas_nB_CpG_only` WGBS
reference in this pipeline was built from (traced via the source archive's own provenance note
pointing at this run under a since-superseded scratch path).

## Cross-pipeline profiling

Stage 5 is meant to be run against [`../flagship-mnase-5hmc-diffbind/`](../flagship-mnase-5hmc-diffbind/)'s
own sonicated-fraction 5hmC bigWig, to ask whether 5hmC is enriched or depleted right at
methylation canyon edges. A DNA methylation bigWig isn't produced by any stage in this pipeline;
the README comment in `05_profile_boundaries.sh` has a one-line `awk`/`bedGraphToBigWig` recipe to
build one from stage 2's `.sym.10xCov.meth` if you want to profile methylation itself alongside
5hmC.

**Not rebuilt (documented as further work, not core scope):** the original archive also
cross-references canyon/UMR boundaries against the flagship's DiffBind fold-classified peak sets
(`posFold`/`negFold`/`NS` per fraction-pair contrast) via `bedtools closest`, in both direction
orderings (`bedtools_overlaps/` vs `bedtools_opp_overlaps/`, i.e. `-a`/`-b` swapped). That's a
real, richer analysis, but reconstructing its exact distance-binning logic (`inside1kb/` subdirs)
would involve more guessing than the core canyon-calling chain above — left as a documented
extension rather than built out.

## Environment

```bash
conda create -n wgbs-canyons -c bioconda -c conda-forge \
  fastqc trim-galore bismark samtools dnmtools bedtools deeptools
```

dnmtools scripts here were checked against dnmtools **1.4.4**'s own CLI usage text (`dnmtools
format`/`uniq`/`counts`/`sym`/`hmr` with no arguments print full usage).
