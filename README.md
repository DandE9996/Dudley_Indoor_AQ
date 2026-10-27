# Indoor environmental quality and retrofit in occupied homes in Dudley, England: analysis code

Analysis code for the article

> **Indoor environmental quality and associated changes after energy retrofit in occupied homes in Dudley, England.**
> J. Qian, E. Pitcher, J. Brooks, H. Langley, S. E. Bartington, H. Li, B. Liu, Y. Shan, Y. Dai, Z. Shi.
> Citation details will be added on publication.

With the required study inputs, the code reproduces the main Tables 1–3 and Figures 1–4, and
Supplementary Tables S1–S6 and Figures S1–S2.
It covers calibrated hourly PM2.5, CO2, temperature, and relative humidity (RH) from 23 occupied homes
(January 2024 to December 2025), and seven retrofit events in four homes. Each event is analyzed on its
own with a difference-in-differences (DiD) model and a controlled interrupted time series (CITS) model,
using never-retrofitted homes monitored over the same dates as comparisons. No pooled retrofit effect is
estimated.

Data processing, descriptive analyses, tables, and figures are Jupyter notebooks (Python). The
statistical models are R scripts that use base R only. **This repository contains code only; no study
data are included** (see [Data](#data)).

## Repository layout

Folders follow the analysis stages. The number at the start of each file name is its position in the
run order.

```
├── 00_shared_functions.ipynb          shared paths, study constants and data loaders (loaded by every notebook)
├── analysis_settings.json             fixed analysis specification (thresholds, windows, models)
│
├── 1_data_conversion/
│   └── 01_convert_hourly_data.ipynb   monitor export (wide) -> one row per sensor, variable and hour
├── 2_calibration/
│   └── 02_calibrate_pm25.ipynb        apply the saved low- and high-concentration PM2.5 calibration models
├── 3_data_cleaning_merging/
│   └── 03_prepare_cohort.ipynb        remove calibration periods, screen 34 monitors, one monitor per home
│                                      (23 homes), merge with home IDs, daily means
├── 4_descriptive_analysis/
│   ├── 04_descriptive_statistics.ipynb  seasonal and diurnal summaries, threshold burdens (Table S1)
│   └── 05_questionnaire_tables.ipynb    household and housing characteristics (Tables 1 and S4)
├── 5_did_cits/
│   ├── 06_event_eligibility.ipynb     event and comparison-home eligibility
│   ├── 07_model_functions.R           REML mixed model (home intercept + AR(1) errors) and helpers
│   ├── 08_fit_did_cits_models.R       56 primary DiD/CITS models + 14 PM2.5 comparison-pool models
│   ├── 09_fit_log1p_sensitivity.R     14 PM2.5 models on the log(1 + x) scale
│   ├── 10_run_placebo_tests.R         treated-label placebo checks
│   └── 11_result_tables.ipynb         Tables 2, 3, S2, S3, S5, S6 and the Figure 3 data
├── 6_figures/
│   ├── 12_Fig1_seasonal_distributions.ipynb
│   ├── 13_Fig2_diurnal_profiles.ipynb
│   ├── 14_Fig3_event_trajectories.ipynb
│   ├── 15_Fig4_event_estimates.ipynb
│   ├── 16_FigS1_threshold_burden.ipynb
│   └── 17_FigS2_PM25_comparison_pool.ipynb
│
├── DATA_REQUIREMENTS.md               the private input files and their formats
├── requirements.txt                   Python packages
├── requirements-calibration.txt       optional H2O package for step 02
└── LICENSE                            MIT License
```

## Workflow

| Step | File | What it does | Article content |
|---|---|---|---|
| 01 | `1_data_conversion/01_convert_hourly_data.ipynb` | Reshapes the hourly monitor export into long format on a complete UTC hourly grid. Missing hours stay empty. | Methods 2.2 |
| 02 | `2_calibration/02_calibrate_pm25.ipynb` | Applies the two saved H2O PM2.5 calibration models (raw PM2.5 below or above 50 µg/m³). Optional; off by default. | Methods 2.2 |
| 03 | `3_data_cleaning_merging/03_prepare_cohort.ipynb` | Removes 31 co-location calibration periods and one relocated monitor, screens 34 monitors (coverage > 35%, CO2/PM2.5 stagnation < 3%), keeps one monitor per residence (23 homes), and forms local-day means from at least 18 valid hours. | Methods 2.2, SI Methods S1 |
| 04 | `4_descriptive_analysis/04_descriptive_statistics.ipynb` | Seasonal distributions, equal-home diurnal profiles and amplitudes, threshold burdens. | Results 3.2, Table S1 |
| 05 | `4_descriptive_analysis/05_questionnaire_tables.ipynb` | Links questionnaires to homes through the sensor register and summarizes them. | Tables 1 and S4 |
| 06 | `5_did_cits/06_event_eligibility.ipynb` | Applies the coverage and PM2.5 baseline rules to every home and event; documents the Home 4 and Home 17 exclusions. | Methods 2.4.1, SI Methods S2 |
| 07 | `5_did_cits/07_model_functions.R` | Model engine and helpers used by steps 08–10 (sourced, not run). | SI Methods S2 |
| 08 | `5_did_cits/08_fit_did_cits_models.R` | Fits the DiD and CITS mixed models for every event and outcome. | Table 3, Fig. 4, Tables S2–S3, Fig. S2 |
| 09 | `5_did_cits/09_fit_log1p_sensitivity.R` | Refits the PM2.5 models with each daily mean x replaced by log(1 + x). | Table S3, SI Methods S3 |
| 10 | `5_did_cits/10_run_placebo_tests.R` | Reassigns the treated label to every participating home and refits. | Table S6, Results 3.3 |
| 11 | `5_did_cits/11_result_tables.ipynb` | Assembles the model result tables and the Figure 3 trajectories. | Tables 2, 3, S2, S3, S5, S6 |
| 12–17 | `6_figures/*.ipynb` | Draws Figures 1–4, S1 and S2 (PNG at 900 dpi, plus PDF or SVG). | Figs. 1–4, S1, S2 |

The PM2.5 calibration models applied in step 02 were developed and evaluated in a companion study
(Qian et al., *Atmospheric Measurement Techniques* 19, 603–615, 2026,
<https://doi.org/10.5194/amt-19-603-2026>). Their training code is archived at
<https://doi.org/10.5281/zenodo.18305879>.

## Requirements

* **Python** 3.12 or later (release validation used Python 3.12) and the packages
  in `requirements.txt`:
  `pip install -r requirements.txt`
* **R** 4.3 or later, base packages only (the article used R 4.6.0). `Rscript` runs steps 08–10.
* Optional, step 02 only: `requirements-calibration.txt` (H2O) and Java.

The figure notebooks use Liberation Sans (Arial is used if it is not installed; the two fonts have the
same character widths). A folder containing `LiberationSans-*.ttf` can be supplied with `INDOOR_AQ_FONT_DIR`.

## Data

The monitoring and questionnaire data could identify participating households and are not public.
De-identified daily data are available from the corresponding author on reasonable request, subject to
the conditions of the ethics approval. Place the files you receive in a folder named `data/` in the
repository (or point `INDOOR_AQ_DATA` to another folder). `DATA_REQUIREMENTS.md` lists every file, its
format, and the step that needs it. Depending on the files available, the workflow can start at step 01
(monitor export), step 03 (calibrated hourly data), or step 04 (prepared hourly and daily data).

The study metadata (home ID to sensor ID crosswalk, retrofit event dates, monitor register, and
co-location calibration periods) are household-level information and are kept with the private data in
`data/` as well. The repository holds only code, documentation, and the fixed analysis settings
(`analysis_settings.json`). `.gitignore` excludes `data/`, `results/`, and data file types.

## How to run

Every step writes to `results/` (or `INDOOR_AQ_OUTPUT`) with the subfolders `data/`, `tables/`,
`figures/` and `logs/`, and reads the outputs of earlier steps from there. Supplied files in `data/` are
never modified. Run the steps in numerical order.

**In Jupyter or VS Code.** Start Jupyter anywhere inside the repository, open each notebook in turn, select
the Python environment with the required packages, and run all cells. Each notebook finds the repository
by searching its own folder and the folders above it. Run the R scripts from a terminal in the repository
folder (`Rscript 5_did_cits/08_fit_did_cits_models.R`, then 09 and 10), or `source()` them in RStudio with
the repository as the working directory.

**From the command line** (in the repository folder; the same commands work in PowerShell and bash):

```
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 1_data_conversion/01_convert_hourly_data.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 2_calibration/02_calibrate_pm25.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 3_data_cleaning_merging/03_prepare_cohort.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 4_descriptive_analysis/04_descriptive_statistics.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 4_descriptive_analysis/05_questionnaire_tables.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 5_did_cits/06_event_eligibility.ipynb
Rscript 5_did_cits/08_fit_did_cits_models.R
Rscript 5_did_cits/09_fit_log1p_sensitivity.R
Rscript 5_did_cits/10_run_placebo_tests.R
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 5_did_cits/11_result_tables.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 6_figures/12_Fig1_seasonal_distributions.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 6_figures/13_Fig2_diurnal_profiles.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 6_figures/14_Fig3_event_trajectories.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 6_figures/15_Fig4_event_estimates.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 6_figures/16_FigS1_threshold_burden.ipynb
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 --output-dir results/notebooks 6_figures/17_FigS2_PM25_comparison_pool.ipynb
```

Start lower in the list when the corresponding prepared files are supplied in `data/` (see
`DATA_REQUIREMENTS.md`). Steps 01–02 are skipped when calibrated hourly data are supplied; step 03 is
skipped when `hourly_primary_23.csv.gz` and `daily_home.csv` are supplied. Runtime depends on hardware, data size, and figure resolution. The command-line examples allow
15 minutes per notebook cell; increase this limit if needed on a slower machine.

### Optional settings (environment variables)

| Variable | Purpose | Default |
|---|---|---|
| `INDOOR_AQ_DATA` | Folder with the supplied data | `data/` |
| `INDOOR_AQ_OUTPUT` | Results folder | `results/` |
| `INDOOR_AQ_CODE` | Repository folder, when Jupyter or R starts outside it | found automatically |
| `INDOOR_AQ_FIGURE_DPI` | Resolution of the PNG figures | 900 |
| `INDOOR_AQ_FONT_DIR` | Folder containing `LiberationSans-*.ttf` | system fonts |
| `INDOOR_AQ_RUN_CALIBRATION` | Set to `1` to score PM2.5 with the H2O models in step 02 | off |
| `INDOOR_AQ_CALIBRATION_MODELS` | Folder with `modCalib_lowlvl` and `modCalib_highlvl` | `data/models/` |
| `INDOOR_AQ_RAW_INPUT` | Uncalibrated long hourly file for step 02 | step 01 output |
| `INDOOR_AQ_CALIBRATED_INPUT` | Calibrated long hourly file for step 03 | step 02 output, otherwise `data/` |

Relative paths in these variables are resolved against the repository folder. The R scripts also accept
the input and output folders as command-line arguments (relative to the current folder); see the header
of each script.

## Analysis specification

The fixed settings are recorded in `analysis_settings.json`.

* Unit of analysis: one home and one retrofit event. Seven events, four outcomes (temperature, RH, PM2.5,
  CO2), and two models give 56 primary models.
* Window: 90 days before works begin and days 1–30 after completion. Construction days are excluded from
  treated and comparison series.
* A daily mean needs at least 18 valid hours on a local (Europe/London) calendar day. Each home needs at
  least 63 valid baseline days and 21 valid follow-up days, and each event needs at least five comparison
  homes for every outcome.
* Comparison homes are never-retrofitted homes. Homes 4 and 17 are excluded from the comparison pool
  because both had recorded works.
* PM2.5 only: the treated home and each comparison home need a baseline mean of at most 100 µg/m³
  for that event. A sensitivity analysis repeats the PM2.5 models with the expanded coverage-eligible
  comparison pool. This threshold was fixed after exploratory analyses and was not prospectively
  prespecified.
* DiD: `y ~ post * treated`. CITS: `y ~ post * treated + time * treated + posttime * treated`, with time
  in 30-day units and `posttime` centered at follow-up day 15.5. A heating-season term (November–March)
  is added when both seasons occur in the window.
* Both models use a home random intercept and AR(1) errors over the actual calendar-day gaps, fitted by
  restricted maximum likelihood from three starting values (the best converged fit is kept). Approximate
  95% intervals use t critical values with G − 2 degrees of freedom, where G is the number of homes.
* The log sensitivity changes only the response, to log(1 + x). Homes, dates, design, and eligibility are
  those of the primary analysis.
* Placebo check: the treated label is moved to each participating home in turn, and the observed absolute
  estimate is ranked among the refits.

## Reproducibility check

Release validation on 27 October 2026 used Windows, Python 3.12.4, NumPy 2.4.2, pandas 2.2.2,
Matplotlib 3.10.8, openpyxl 3.1.4, and R 4.6.0. All 14 notebooks ran from top to bottom in fresh kernels,
and scripts 08–10 refitted the 70 primary/expanded-pool models, 14 log1p models, and 57 placebo
distributions. The six figure notebooks generated PNGs at 900 dpi.

This run started from the monitor export for step 01 and used the supplied, previously calibrated
hourly data for step 03. Step 02 ran with H2O scoring disabled; applying the saved H2O models was not
retested in this release check. Private input files were unchanged, and all generated results were
written to a separate validation directory outside the public code folder.

Compared with the saved study outputs, all 65,268 daily records matched in keys, valid-hour counts,
and missingness; the maximum absolute difference in daily means was 5.7 × 10⁻¹⁴. All 70 DiD/CITS
models converged and retained the same comparison counts. Maximum differences were 1.1 × 10⁻⁵ in
estimates and 1.6 × 10⁻⁴ in interval limits; every estimate and interval limit agreed at two decimal
places. All 14 log1p models converged, with differences below 5 × 10⁻¹³. All 56 primary placebo ranks
and group sizes matched the saved reference, and none of the refits in the 57 placebo distributions
failed. The built-in independent dense-GLS checks also passed.

These checks used authorized private study inputs. A public checkout provides the analysis code and
instructions, but cannot reproduce the results until the required inputs are supplied.

## License

The code and accompanying documentation are released under the [MIT License](LICENSE).
Study data and saved calibration models are not included in this repository; the code license does
not grant access to those files.
