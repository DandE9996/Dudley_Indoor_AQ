# Data requirements

No study data are part of this repository, because they could identify participating households.
This includes the study metadata below (home crosswalk, event dates, monitor register, calibration
periods). De-identified daily data are available from the corresponding author on
reasonable request, subject to the conditions of the ethics approval.

Place supplied files in `data/` in the repository folder, or set `INDOOR_AQ_DATA` to another folder.
Files produced by an earlier step in `results/data/` take precedence over files of the same name in
`data/`. Nothing in `data/` is modified.

## Which files each step needs

For a complete run from step 03, supply the study metadata files `home_id_crosswalk_23.csv`,
`events.csv`, `registered_sensors.csv`, and `calibration_intervals.csv` in `data/`.
The measurement files depend on where the workflow starts:

| Start the workflow at | Files needed in `data/` |
|---|---|
| Step 01 (monitor export) | `indoor_all.csv`, plus `indoor_all_long_hourly_calibrated.csv` unless step 02 scores the data |
| Step 03 (calibrated hourly data) | `indoor_all_long_hourly_calibrated.csv` |
| Step 04 (prepared data) | `hourly_primary_23.csv.gz` and `daily_home.csv` |
| Step 05 (questionnaire tables) | `Brockmoor_home_info.xlsx`, `Survey.xlsx`, `Dudley-NZN-Sensor-Retrofit-Data - UoB.xlsx` |

Step 05 is independent of the other analysis stages and can be skipped when the questionnaires are
not available. Steps 00 and 07 only define or inspect the shared configuration; they can be run without
study measurements. Step 02 can also be inspected with scoring disabled.

For a partial run, the actual dependencies are:

| Steps | Required inputs |
|---|---|
| 01 | `indoor_all.csv`; no study metadata CSVs are read |
| 02, scoring enabled | Uncalibrated long hourly data and both saved H2O models |
| 03 | Calibrated long hourly data and all four study metadata CSVs |
| 04 and 06 | Prepared hourly and daily data plus the home crosswalk; step 06 also reads `events.csv` |
| 05 | The three questionnaire/register workbooks and the home crosswalk |
| 08 | `daily_home.csv`, `events.csv`, and the home crosswalk |
| 09–10 | `data/event_models.rds` produced by step 08 in the selected results folder |
| 11 | Model, observed-mean, and analysis-row outputs from step 08; log-sensitivity output from step 09; placebo ranks from step 10; and `events.csv` in the private input folder |
| 12–13 | Prepared hourly data and the home crosswalk |
| 14 | `trajectory_plot_data.csv`, `event_results.csv`, and the copied `data/events.csv` from this run |
| 15 and 17 | `event_results.csv` and the copied `data/events.csv` from this run |
| 16 | Prepared hourly and daily data plus the home crosswalk |

Prepared measurements alone do not replace the model outputs required by steps 09–11, 14–15, and 17.
Run the corresponding earlier stages before executing those steps. Generated tables and model files
are not distributed with this repository.

## File formats

### Study metadata

| File | Columns used | Content |
|---|---|---|
| `home_id_crosswalk_23.csv` | `home_id`, `home`, `source_home`, `sensor_id`, `is_treated` | Home 1 to Home 23 and the retained monitor of each home |
| `events.csv` | `event`, `home`, `start`, `end` | The seven analyzed retrofit events (`H05_Loft`, `H05_DMEV`, `H10_Loft`, `H10_DMEV`, `H14_CWI`, `H14_Loft`, `H16_Bundle`), first and last day of works |
| `registered_sensors.csv` | `sensor_id` | The 34 registered monitors screened in step 03 |
| `calibration_intervals.csv` | `sensor`, `start`, `end`, `prior_audit_hours` | The 31 co-location calibration periods (UTC, inclusive) removed in step 03, with their audited hour counts |

Other columns are ignored.

### `indoor_all.csv` (step 01)

Hourly export of the AirGradient ONE monitors, one row per monitor and hour.

| Column | Content |
|---|---|
| `Location Name` | Monitor name ending in its numeric sensor ID (for example `UOB_009` is sensor 9) |
| `UTC Date/Time` | Start of the hour in UTC (ISO 8601) |
| `PM2.5 (μg/m³) corrected`, `CO2 (ppm) corrected`, `Temperature (°C) corrected`, `Humidity (%) corrected` | Hourly values; empty when missing |

Other columns are ignored. Columns are located by their prefixes (`PM2.5`, `CO2`, `Temperature`, `Humidity`).

### `indoor_all_long_hourly_calibrated.csv` (step 03)

Long format, one row per sensor, variable, and hour, on a complete hourly grid. Step 01 writes the same
layout without calibration (`indoor_all_long_hourly.csv`), and step 02 writes this file when scoring is
enabled.

| Column | Content |
|---|---|
| `sensor_id` | Numeric sensor ID |
| `variable` | `co2`, `pm25`, `rh`, or `temp` |
| `time` | Hour in Europe/London time with UTC offset (for example `2024-06-01 13:00:00+01:00`) |
| `value` | Measurement (PM2.5 calibrated); empty when missing |

Step 03 checks the input against the study cohort: 34 registered monitors, 31 co-location calibration
periods with their audited hour counts, 25 monitors passing screening, and the 23 homes in
`home_id_crosswalk_23.csv`. A different input stops with an assertion message.

### `models/modCalib_lowlvl` and `models/modCalib_highlvl` (step 02, optional)

The saved H2O models for raw PM2.5 below and at or above 50 µg/m³. Set `INDOOR_AQ_RUN_CALIBRATION=1`
to score, and `INDOOR_AQ_CALIBRATION_MODELS` if the models are stored elsewhere. Scoring needs
`requirements-calibration.txt` and Java. The models and their training code are described in
Qian et al. (2026), <https://doi.org/10.5194/amt-19-603-2026>, with code at
<https://doi.org/10.5281/zenodo.18305879>.

### `hourly_primary_23.csv.gz` and `daily_home.csv` (prepared data, written by step 03)

`hourly_primary_23.csv.gz`: hourly values of the 23 retained monitors.

| Column | Content |
|---|---|
| `sensor` | Numeric sensor ID (homes are resolved through `home_id_crosswalk_23.csv`) |
| `home` | Home label (informational; the crosswalk is authoritative) |
| `time` | Hour in UTC (ISO 8601) |
| `variable` | `co2`, `pm25`, `rh`, or `temp` |
| `value` | Calibrated measurement; empty when missing |

`daily_home.csv`: one row per home, variable, and local (Europe/London) calendar day.

| Column | Content |
|---|---|
| `sensor_id` | Numeric sensor ID |
| `home` | Home label, `Home 1` to `Home 23`, as in the crosswalk |
| `variable` | `co2`, `pm25`, `rh`, or `temp` |
| `date` | Local date (`YYYY-MM-DD`) |
| `n_hours` | Number of valid hourly values |
| `y` | Daily mean; empty when fewer than 18 valid hours |

### Questionnaire workbooks (step 05)

| File | Content | Layout used by the code |
|---|---|---|
| `Brockmoor_home_info.xlsx` | Sensor register linking each monitor to a residence | First worksheet; address ID in column A (merged cells allowed), sensor ID in column B |
| `Survey.xlsx` | Dwelling and occupancy assessment (Table 1) | First worksheet; header in row 1, including `Address ID` |
| `Dudley-NZN-Sensor-Retrofit-Data - UoB.xlsx` | Project housing survey (Table S4) | Worksheet `Survey Responses`; header in row 2, including `Address ID` |

Records are linked through the address ID in the sensor register, never by row order. The question
columns read from each survey are listed in `05_questionnaire_tables.ipynb`.
