# Portfolio Pipelines

Cleaned-up, runnable versions of genomics analysis pipelines originally developed for wet/dry
lab research, repackaged here for demonstration purposes. Each pipeline is organized as a
sequence of self-contained, single-purpose stage scripts with explicit inputs/outputs/params
declared in a header comment, so each stage can be run standalone, dropped into a SLURM array,
or later wrapped in Nextflow/Snakemake without a rewrite.

Demo inputs use public GEO/SRA accessions wherever available (see each pipeline's README for
the specific accession(s) and how to pull the data). Where public data can't stand in for the
original inputs, that is documented rather than faked.

## Pipelines

| Pipeline | Stages | Description |
|---|---|---|
| [`flagship-mnase-5hmc-diffbind/`](flagship-mnase-5hmc-diffbind/) | 5 | MNase-fractionated chromatin profiling (mono/submono/sonicated) + DiffBind differential-occupancy classification |
| [`chromhmm-20state/`](chromhmm-20state/) | 4 | ChromHMM 20-state chromatin annotation from a 46-mark ChIP-seq panel |
| [`wgbs-dnmtools-canyons/`](wgbs-dnmtools-canyons/) | 5 | WGBS/DNMTools 5hmC–5mC methylation canyon and HMR (hypomethylated region) calling |

All three are fully built and tested against their real tools (not just syntax-checked) — see
each pipeline's README for what "tested" meant given none of the underlying tools have unit-test
suites of their own: synthetic smoke tests through the actual installed binaries, not full runs
against the demo data.

## Running one end-to-end

Every stage script is self-contained: run `./scripts/NN_stage_name.sh` with no arguments to see
its usage, or read the header comment for the full INPUTS/OUTPUTS/PARAMS contract. There's no
top-level runner by design — string the stages together yourself (a shell loop, a SLURM array, or
eventually a Nextflow/Snakemake wrapper) so each stage stays inspectable on its own. The general
shape, common to all three pipelines:

1. Pull the pipeline's demo FASTQ(s) from the public accession in its README's "Demo data" section.
2. Run the align/QC stage(s) per sample.
3. Run the peak-calling / methylation-calling stage(s) per sample.
4. Run the final aggregation stage (DiffBind / ChromHMM LearnModel / canyon derivation) across all
   samples at once.

See "What a full-scale run actually needs" below before running any pipeline against its full
demo dataset rather than a one- or two-sample smoke test.

## What a full-scale run actually needs

These pipelines were built and verified with synthetic, single-chromosome smoke tests — enough
to prove each script's logic is correct, not to characterize real runtime/memory/storage at full
scale. Rough figures below are derived from each stage's actual tool/algorithm and the demo
dataset's real size (read counts, mark counts, etc.), not measured end-to-end on this hardware.

### flagship-mnase-5hmc-diffbind

- **Data volume:** 12 public runs (`SRP684244`), ~20-40M paired reads each, mouse (mm10).
- **Per-sample compute:** stage 01 (trim+align+dedup) is the heavy step — Bowtie2 `--very-sensitive`
  on ~30M PE reads is typically 1-3 hours on 8 cores; stage 02 (T4 spike-in) is much faster (small
  genome); stage 03 (MACS2) is minutes. Run stages 01-03 per sample — 12 samples, each independent,
  so trivially parallel across a SLURM array.
- **Storage:** raw+trimmed FASTQs + BAMs per sample run ~5-15GB; ×12 samples is roughly 100-150GB
  if you keep every intermediate. Delete `.bam` (unsorted) and `.filt.bam`/`.s.filt.bam`
  intermediates once `.s.rms.filt.bam` exists, per sample, to cut that substantially.
  Stage 04 (DiffBind, in-memory over 12 BAMs' worth of counts) needs a few GB of RAM, not much
  disk beyond its bed/tsv outputs.
- **What's missing to run for real:** a `mono`-fraction dataset — the public accession only has
  submono+sonicated, so `comp1`/`comp2` (any contrast involving `mono`) will just log a skip
  message and produce no output until you supply that fraction's BAMs/peaks yourself.

### chromhmm-20state

- **Data volume:** 46 marks (45 public via `SRP075985` + 1 from the flagship's own output),
  single-end, ~25-30M reads each, plus one shared input — this is the largest data volume of the
  three by mark count, even though each individual run is smaller than the flagship's PE data.
- **Per-sample compute:** stage 01 (SE trim+align+dedup) is ~30-60 min/mark on 8 cores; stage 02
  (MACS2) is minutes. 46 independent per-mark jobs — a SLURM array of 46 is the natural shape.
  **Stage 04 (LearnModel) is the real bottleneck**: fitting a 20-state HMM over 46 marks
  genome-wide (mm10, ~2.7Gb at 200bp bins ≈ 13M bins) took hours in the original run and is
  CPU + memory-bound on iteration count, not I/O; give it as many `-p` processors as you can and
  expect a multi-hour job even with everything else already done.
- **Storage:** ~46 × 3-8GB (FASTQ+BAM) ≈ 150-350GB if kept uncompressed/undeleted; binarized
  genome-wide matrices (stage 03) add a few more GB; LearnModel's own output (POSTERIOR/ per
  chromosome, segmentation, plots) is comparatively small (<1GB).
- **What's missing to run for real:** nothing — this is the one pipeline with full public
  coverage of its mark panel (see its README), aside from needing the flagship's own sonicated
  5hmC peaks for the `Son_5hmC` mark first.

### wgbs-dnmtools-canyons

- **Data volume:** one run (`SRR1003257`), single-end WGBS, ~1.27B reads — by far the largest
  single input of the three pipelines (whole-genome bisulfite sequencing at real depth is
  intrinsically large).
- **Per-sample compute:** Bismark alignment of ~1.27B SE reads is the dominant cost — realistically
  many hours to over a day even with `--parallel`, since Bismark's own aligner (Bowtie2 under the
  hood, 3-letter-converted genome) is slow relative to a normal DNA aligner; this is the stage
  most worth testing on a chromosome subset first (align to a single chromosome's Bismark index) before
  committing to a full run. Stages 02-04 (dnmtools chain, HMR calling, canyon derivation) are
  comparatively fast — minutes given the aligned BAM.
- **Storage:** the raw FASTQ is ~55GB compressed (verified against ENA's own file report for this
  run, not estimated) — trimmed FASTQ, the Bismark BAM, and its dedup/sorted copies each add
  comparable amounts transiently. Budget 300-500GB of scratch space for this one run if you keep
  every intermediate; delete the pre-dedup BAM once `.deduplicated.bam` exists.
- **What's missing to run for real:** nothing dataset-wise (the public run covers the whole
  pipeline); a combined single-FASTA mm10 reference for `dnmtools counts -c` (per-chromosome
  FASTAs need concatenating first — see stage 02's header comment).

### Common to all three

- **Accounts/access:** an SRA/ENA-reachable network path (`prefetch`/`fasterq-dump` or direct FTP)
  to pull the FASTQs — no login required, these are all public accessions.
- **Reference genomes:** mm10 Bowtie2 index (flagship, chromhmm), mm10 Bismark-prepared genome
  folder (wgbs), and for the T4 spike-in stage, a Bowtie2 index built from any T4 phage FASTA
  (e.g. RefSeq `NC_000866.4`) — none of these are bundled in this repo.
- **A scheduler helps but isn't required:** every script takes explicit flags rather than assuming
  SLURM, so the same commands work as a bash loop on a workstation, just slower wall-clock than a
  cluster array running samples in parallel.

## Provenance

These pipelines are adapted from analyses supporting:

> Hsu, Chen, Lay, Ma, Delatte, Lio. "The majority of 5hmC is nucleosome-associated in naïve B
> cells." *BMC Genomics* 2026. DOI: [10.1186/s12864-026-13032-y](https://doi.org/10.1186/s12864-026-13032-y)

Scripts are cleaned up and generalized from the original SLURM/HPC batch scripts; where an
original driver script no longer exists (only downstream outputs survived), the README for that
stage says so explicitly and documents what was reconstructed and from what evidence.

## License

MIT — see [LICENSE](LICENSE).
