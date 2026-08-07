# grothlab/crepas: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.1.0dev](https://github.com/grothlab/crepas/releases/tag/1.1.0dev)]

Development version of grothlab/crepas.

### `Added`

- Full `meta.yml` documentation, `nf-test` tests, and topic-based version reporting for the local modules that previously lacked them.
- Comprehensive stubs for all local modules, enabling end-to-end dry-runs of the pipeline (`-stub`) ([#61](https://github.com/grothlab/crepas/issues/61)).
- [`minibwa`](https://github.com/lh3/minibwa) as an additional read-alignment option.
- GEFION HPC configuration profile.

### `Changed`

- Replaced generic `ubuntu:22.04` module containers with purpose-built Wave / BioContainers images.
- Migrated version reporting from per-module `versions.yml` files to the `versions` topic channel.
- Replaced several local modules with their nf-core equivalents (`gtf2bed` now uses ea-utils, plus `SAMTOOLS_REHEADER` and `BEDTOOLS_GENOMECOV`).
- Updated the Chromap nf-core module and reworked `fasta`/`fai` input-channel handling.
- Updated the pipeline to the nf-core tools 4.0.2 template.

### `Fixed`

- flT3 orphan removal is no longer skipped, so TE counting can run on pre-blacklist-filtered BAMs.
- Chromap segmentation-fault and memory issues.
- TEtranscripts and TElocal container definitions.
- `test_cutandrun` was mapped against the wrong genome.
- Missing `versions` workflow output parameter.

## [[1.0.0](https://github.com/grothlab/crepas/releases/tag/1.0.0)] - Mercurian Cinnabar - 2026-06-21

Initial release of grothlab/crepas.
