# Distribution notes

Source snapshot: local revision, 2026-09-10. The PDF and manuscript source are the revised review candidate. Computational caches, per-replicate RDS files, object/shared-library binaries, and local application data were excluded. Aggregate CSVs, source code, tests, reports and figure PDFs are included. The application reads the de-identified merged dataset already in this repository.

Only path configuration was adapted in the distributed R scripts, using LRC_REVISION_ROOT (default revision from repository root):
- 05-writing/make_tables.R
- 05-writing/make_fig_sc4.R
- 02-validation/sc4_boundary/run_sc4_boundary.R
- 04-application/code/prep_data.R
- 04-application/code/make_figures.R
- 04-application/code/build_table.R
- 04-application/code/run_lrcbart.R
- 04-application/code/run_reference.R

The obsolete local `chain_BC.sh` scheduler was omitted because it waits on a machine-specific process ID; the individual full-study R drivers are included. Existing source whitespace and editorial diff records were preserved.
