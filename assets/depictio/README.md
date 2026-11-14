# Depictio template for grothlab/crepas

A [Depictio](https://depictio.github.io/depictio-docs/stable/) template: it turns a crepas
`--outdir` into a Depictio project with an interactive dashboard. Each successful run copies it,
with the files it reads, to `<outdir>/depictio/`; see "Interactive report with Depictio" in
`docs/usage.md` for viewing that on your own computer. By hand, from a full output directory:

```bash
depictio ingest <outdir> --template <path to crepas>/assets/depictio
```

| File                   | What it is                                                                                 |
| ---------------------- | ------------------------------------------------------------------------------------------ |
| `template.yaml`        | What to ingest: the data collections (tables scanned from the output directory) and links  |
| `dashboards/base.yaml` | What to show: the Overview, Reads & filtering, Peaks and MultiQC tabs                      |

## How it reads the results

Depictio runs Python recipes only from its own package, so this template has none. Every table is
a scan of a file the pipeline writes, and the parts of its path that matter (aligner, experiment
type, peak caller, library) become columns through scan wildcards. The one file written for the
template is the library hub, `pipeline_info/libraries.tsv`; every per-library table and the
MultiQC panels are linked to it on the library id, which is what lets the library filters narrow
every tab.

All collections except the hub are optional, and a dashboard tile whose collection or MultiQC
module is missing from a run is dropped at import, so one template serves every experiment type,
aligner and set of peak callers.

At the end of each successful run, `gatherDepictioReport` (in
`subworkflows/local/utils_grothlab_crepas_pipeline/main.nf`) copies this folder to
`<outdir>/depictio/template/` and every file the template reads to `<outdir>/depictio/data/`. It
works out which files from this `template.yaml` (scan patterns, single-file paths and provenance
globs), so a new or changed collection needs no change there.

## Changing it

- Tags must agree between the two files: `workflow_tag: crepas`, and each `data_collection_tag`.
- Bump `template.version` with each change: major when a collection or a column changes (projects
  must be ingested again), minor for dashboard-only changes, patch for text.
- Check a dashboard edit with `depictio dashboard validate dashboards/base.yaml --offline`, then
  ingest a real output directory into a throwaway server (`depictio local up`, on your own
  machine or inside a job, never on a shared login node) and look at every tab.
- Code-mode figures get a Polars `df` and `px`/`pl`, take no imports, and must build their data as
  `df_modified` before the `fig = px...` line.
- If a crepas output file or directory is renamed, update the matching `pattern` here.
