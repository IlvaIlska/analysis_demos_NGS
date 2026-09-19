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

| Pipeline | Status | Description |
|---|---|---|
| [`flagship-mnase-5hmc-diffbind/`](flagship-mnase-5hmc-diffbind/) | Built out | MNase-fractionated chromatin profiling (mono/submono/sonicated) + DiffBind differential-occupancy classification |
| [`chromhmm-20state/`](chromhmm-20state/) | Scaffolded | ChromHMM 20-state chromatin annotation from a ~35-mark ChIP/CUT&Tag panel |
| [`wgbs-dnmtools-canyons/`](wgbs-dnmtools-canyons/) | Scaffolded | WGBS/DNMTools 5hmC–5mC methylation canyon and HMR (hypomethylated region) calling |

## Provenance

These pipelines are adapted from analyses supporting:

> Hsu, Chen, Lay, Ma, Delatte, Lio. "The majority of 5hmC is nucleosome-associated in naïve B
> cells." *BMC Genomics* 2026. DOI: [10.1186/s12864-026-13032-y](https://doi.org/10.1186/s12864-026-13032-y)

Scripts are cleaned up and generalized from the original SLURM/HPC batch scripts; where an
original driver script no longer exists (only downstream outputs survived), the README for that
stage says so explicitly and documents what was reconstructed and from what evidence.

## License

MIT — see [LICENSE](LICENSE).
