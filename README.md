# Replication package

**Children's Education and Parental Longevity. Causal Evidence From Two English School-Leaving Age Reforms**

This package reproduces every table, figure and in-text statistic computed from data in the paper and its supplement. Table 1 (literature summary) and Table A1 (ethical approval) are compiled by hand.

## Requirements

- R 4.5.3
- R packages: tidyverse 2.0.0, haven 2.5.5, survival 3.8.6, survminer 0.5.2, patchwork 1.3.2, car 3.1.5, rddensity 3.0, rdrobust 4.0.0, mice 3.19.0
- About one hour on a standard desktop computer

## Data

Place the files as follows (paths relative to the package root).

| Source | Used for | Location |
|---|---|---|
| English Longitudinal Study of Ageing, waves 1-3 (SN 5050) and waves 4-11 (SN 8346), UK Data Service, special licence, Stata files | main text, Supplements B-E | `data/elsa sn 5050/`, `data/elsa sn 8346/` |
| ELSA waves 1-8, NatCen Social Research restricted-access version, SPSS files `Wave1_final.sav` ... `Wave8_final.sav` | Supplements D and E | `data/elsa from natcen/` |
| Health Survey for England 1991-2013, UK Data Service, Stata files | Figure 1 | `data/hse/` |

Scripts `00`, `02` and `03` state the exact file names and sub-folders they read.

## Running the code

From the package root:

```bash
Rscript run_all.R
```

`run_all.R` runs the scripts in `code/` in order, each in its own R process, and writes all R output to `logs/replication_log.txt`. A banner before each step names the table or figure it produces. Tables are written to `Outputs/` (CSV) and figures to `figures/` (PDF). Intermediate data files are written to `data/tempdata/`.

A single script can be run with `Rscript code/<script>` from the package root once scripts 00-02 have run.

## Contents

| Script | Produces | Output files |
|---|---|---|
| `00_prepare_ukds_data.R` | Base data set from SN 5050 / SN 8346 | `data/tempdata/ELSA_base_dataset.rds` |
| `01_prepare_ukds_analysis_file.R` | Analysis file, UK Data Service | `data/tempdata/ELSA_sn8346_bmonth-year.rds` |
| `02_prepare_natcen_data.R` | Analysis file, NatCen | `data/tempdata/ELSA_bmonth-year.rds` |
| `03_figure1.R` | Figure 1 | `figures/figure1.pdf` |
| `04_table2_figures2to4.R` | Table 2, Panels A and C; Table A4; Table A5, Panels A and C; Table A6, UK Data Service columns; Table A8, IV column; Figures 2, 3, 4 | `Outputs/table2.csv`, `first_stage.csv`, `figure4_estimates.csv`; `figures/figure2.pdf`, `figure3.pdf`, `figure4.pdf` |
| `05_table2_panelB.R` | Table 2, Panel B; Table A5, Panel B | `Outputs/table2_panelB.csv` |
| `06_tablesA2_A3.R` | Tables A2, A3 | `Outputs/tableA2.csv`, `tableA3.csv` |
| `07_figuresA1toA4.R` | Figures A1-A4 (Supplement C.1-C.4) | `figures/figureA1.pdf` ... `figureA4.pdf`; `Outputs/schoenfeld.csv`, `balance.csv`, `rddensity_1947.txt`, `rddensity_1972.txt`, `figureA4_estimates.csv` |
| `08_proportional_hazards_C1.R` | Supplement C.1 text statistics | `Outputs/ph_timevarying.csv`, `ph_strata.csv`, `ph_bands.csv` |
| `09_missing_data_weights_figureA5.R` | Figure A5 (Supplement C.5) | `figures/figureA5.pdf`; `Outputs/mice.csv`, `weights.csv`, `figureA5_estimates.csv` |
| `10_tableA6_natcen.R` | Table A6, NatCen columns | `Outputs/tableA6_natcen.csv` |
| `11_sample_overlap_D1.R` | Supplement D.1 text statistics | log only |
| `12_tableA7.R` | Table A7 | `Outputs/tableA7.csv` |
| `13_tableA8.R` | Table A8, predicted-hazard RDD column | `Outputs/tableA8.csv` |
| `helpers.R` | Variable construction, first stage and bootstrap functions shared by scripts 04-13 | |

Tables in the manuscript are typeset by hand from the CSV files.

## Reproducibility notes

- Seeds: `set.seed(321)` in scripts 04 and 10; `set.seed(20261008)` in script 09, where each multiple-imputation run is additionally seeded with `seed = 1`.
- Bootstrap replications: 1,000 (Table 2, Figure 4, Tables A5 and A6); 100 per imputed data set with 20 imputations (Figure A5, Panel A); 500 (Figure A5, Panel B).
- Scripts read their inputs from `data/tempdata/` and must be run in numerical order.
- Warnings printed by R are kept in the log.
