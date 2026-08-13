# grothlab/crepas: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.1.0dev](https://github.com/grothlab/crepas/releases/tag/1.1.0dev)]

Development version of grothlab/crepas.

### `Added`

- Full `meta.yml` documentation, `nf-test` tests, and topic-based version reporting for the local modules that previously lacked them.
- Comprehensive stubs for all local modules, enabling end-to-end dry-runs of the pipeline (`-stub`) ([#61](https://github.com/grothlab/crepas/issues/61)).
- [`minibwa`](https://github.com/lh3/minibwa) as an additional read-alignment option.
- `hybrid_fasta` parameter and improved spike-in genome handling.
- Per-run provenance / RO-Crate ([WRROC](https://www.researchobject.org/workflow-run-crate/)) generation via the [`nf-prov`](https://github.com/nextflow-io/nf-prov) plugin.
- GEFION HPC configuration profile.
- Support for the `eSPAN` experiment type ([#28](https://github.com/grothlab/crepas/issues/28)). Includes a `test_espan` profile and test samplesheet.
- Grouped partition plots, combining all samples sharing an antibody into a single plot, with a `skip_partition_group_plot` parameter to disable them.
- A reduced `*.totals.tsv` summary table emitted alongside the full samtools stats summary, holding only the per-processing-step read counts.
- Experiment-type-specific alignment settings for CUT&RUN, CUT&Tag, TIP-seq and eSPAN, applied to every aligner that exposes the corresponding options (bowtie2, bowtie, HISAT2, Chromap and minimap2).
- `peak_callers`, a comma-separated list of peak callers that are all run in the same pipeline run (e.g. `--peak_callers macs3,genrich`), replacing the single-valued `peak_caller`. The requested caller is carried in each sample's `meta`, so the results of different callers no longer overwrite each other.
- `skip_peak_calling`, to switch peak calling off entirely without changing the `peak_callers` selection.
- `partition_iz_rm_overlap_range`, which sets the window used to discard overlapping initiation zones independently of `partition_plot_range`, so that the plotted range can be widened without discarding more initiation zones.

### `Changed`

- Replaced generic `ubuntu:22.04` module containers with purpose-built Wave / BioContainers images.
- Migrated version reporting from per-module `versions.yml` files to the `versions` topic channel.
- Replaced several local modules with their nf-core equivalents (`gtf2bed` now uses ea-utils, plus `SAMTOOLS_REHEADER` and `BEDTOOLS_GENOMECOV`).
- Updated the Chromap nf-core module and reworked `fasta`/`fai` input-channel handling.
- Removed the pipeline-wide `sort_bam` parameter; alignment sorting, indexing and stats are now always handled by `BAM_SORT_STATS_SAMTOOLS`, with the aligners emitting an unsorted BAM.
- Updated numerous nf-core modules (aligners, deepTools, MACS3, HOMER, bedtools, MultiQC, FastQC, …) to their latest versions.
- Added Apptainer container definitions to the local modules.
- Updated the pipeline to the nf-core tools 4.0.2 template, bumped `nf-schema` to 2.7.3 and `nf-prov` to 1.7.0, and regenerated the parameter documentation from the schema.
- `BAM_SPLIT_BY_STRAND` now takes `exp_type` and `strandedness` as explicit inputs and selects reads with `samtools view --expr` filter expressions, excluding unmapped records from every output.
- Renamed the samtools stats summary outputs to `*.all.tsv` (every column) and `*.totals.tsv` (read counts only), so that the two tables can be matched separately.
- Renamed the multimapping-reads column of the samtools stats summary to `<column>_minus_<column>`, since the filtering step it is derived from removes more than just multimapping reads.
- The GEFION profile now prepends to the container `PATH` instead of replacing it, so images that install elsewhere (e.g. Wave/pixi images) keep working.
- Expanded `CITATIONS.md.
- SEACR is no longer restricted to CUT&RUN, CUT&Tag and TIP-seq samples; it now runs on every sample for which it is requested through `peak_callers`.
- Updated the pipeline metro map.
- The samplesheet examples in `docs/usage.md` now match the example samplesheets shipped in `assets/test-datasets/`, and the DAN System guide (`docs/ku_sund_danhead_crepas_usage.md`) lists every available test profile.
- Renamed `min_reps_consensus` to `consensus_min_replicates_per_sample`, and the corresponding `macs3_merged_expand.py` argument from `--min_replicates` to `--min_replicates_per_sample`.
- Removed some channel dumps written to `<outdir>/.debug/BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2/`.

### `Fixed`

- flT3 orphan removal is no longer skipped, so TE counting can run on pre-blacklist-filtered BAMs.
- Chromap segmentation-fault and memory issues.
- TEtranscripts and TElocal container definitions.
- `test_cutandrun` was mapped against the wrong genome.
- Missing `versions` workflow output parameter.
- Incorrect `fasta`/chromsizes selection when both `hybrid_fasta` and `fasta` were provided.
- `strobealign` ignoring the spike-in genome index.
- Division-by-zero when a sample had no spike-in reads (now guarded and logged).
- Endogenous BAM normalization incorrectly skipped when the exogenous/spike-in BAM had no reads.
- `params.seq_platform` not being applied to the BAM read group.
- Parameter validation guards that tested `params.containsKey()` on parameters declared with a `null` default, so they were unreachable (blacklist, initiation zones / OK-seq RFD file, GTF/GFF and TE index checks). As a result, enabling blacklist filtering without a blacklist silently emptied every downstream channel while the run still reported success.
- CISRPM normalization crashing with `Cannot invoke method div() on null object` when the input control had no spike-in reads: the filter guarded the input control's endogenous total but the normalization factor used its exogenous one.
- Duplicated `--maxins` in the bowtie2 arguments, which silently overrode the CUT&RUN, CUT&Tag and TIP-seq insert-size limit with the default value.
- Bugs in `process_stats_summary.R`.
- An empty `final_samtools_stats_summary` table, and a stale `meta` reference in `STATS_CAT`.
- Add missing `skip_flTbl` param in profiles to avoid.
- `PARTITION_PLOT` failing for samples without an input control.
- Partition and RFD plot labels for sample names containing dots (e.g. `H3.3`), which were truncated at the first dot.
- File name collisions in `GENRICH` when the same input control is shared between several IP samples; treatment and control BAMs are now staged into separate directories.
- Consensus peak `plotProfile` outputs of different experiment types sharing an antibody being written to the same directory; the output path and file prefix now include the experiment type.
- `featureCounts` quantification and DESeq2 QC of the consensus peaks not running.

## [[1.0.0](https://github.com/grothlab/crepas/releases/tag/1.0.0)] - Mercurian Cinnabar - 2026-06-21

Initial release of grothlab/crepas.
