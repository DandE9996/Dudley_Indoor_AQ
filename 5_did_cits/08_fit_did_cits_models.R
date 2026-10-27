# =============================================================================
# 08_fit_did_cits_models.R  (step 08)
#
# Fit the event-specific difference-in-differences (DiD) and controlled
# interrupted time series (CITS) mixed models.
#
# Manuscript : Table 3 and Figure 4 (56 primary models); Tables S2 and S3 and
#              Figure S2 (14 PM2.5 models with the expanded comparison pool)
# Usage      : Rscript 5_did_cits/08_fit_did_cits_models.R [DAILY_DATA_DIR] [OUTPUT_DIR]
#                DAILY_DATA_DIR folder containing daily_home.csv. Default:
#                               results/data (written by step 03), otherwise
#                               data/ (INDOOR_AQ_DATA)
#                OUTPUT_DIR     results folder. Default: results/
#                               (INDOOR_AQ_OUTPUT)
#              events.csv and home_id_crosswalk_23.csv are read from data/
#              (INDOOR_AQ_DATA), with the other private study files.
# Runtime    : under a minute (70 models x 3 optimizer starts)
#
# Study design (also recorded in analysis_settings.json)
# ------------------------------------------------
# * One analysis per treated home x retrofit event x outcome: 7 events x
#   4 outcomes (temp, rh, pm25, co2) x 2 models = 56 primary models.
#   No pooled retrofit effect is estimated.
# * Window: the 90 days before works start and days 1-30 after works end.
#   Construction days are excluded.
# * Comparison homes: never-treated homes, excluding Homes 4 and 17 (which
#   had recorded works). A home is used for an outcome when it has at least
#   63 valid baseline days and 21 valid follow-up days. For PM2.5, the
#   primary analysis also requires a baseline mean <= 100 ug/m3; the SI
#   sensitivity analysis drops that restriction. Every analysis needs at
#   least 5 comparison homes.
# * Fixed effects:
#     DiD  : y ~ post * treated [+ season]
#     CITS : y ~ post * treated + time * treated + posttime * treated [+ season]
#   time = (date - start) / 30 and posttime = post * (date - end - 15.5) / 30.
#   The retrofit effect is the `post:treated` coefficient in both models; in
#   CITS it is the level change at the middle of the 30-day follow-up.
#   `season` (heating season = November-March) is added only when the window
#   contains both seasons.
# * Random effects and inference: home random intercept and AR(1) errors
#   over actual calendar-day gaps, fitted by REML (mixed_fit() in 07_model_functions.R);
#   best of three optimizer starts; approximate t(G - 2) 95% intervals
#   (07_model_functions.R).
#
# Outputs (inside OUTPUT_DIR)
# ---------------------------
#   tables/event_results.csv                 all 70 models
#   tables/primary_event_results.csv         the 56 primary models
#   tables/SI_PM25_original_comparisons.csv  the 14 expanded-pool PM2.5 models
#   tables/home_event_eligibility.csv        inclusion decision for every home
#   tables/home_event_observed_means.csv     observed before/after means
#   tables/optimizer_checks.csv              outcome of every optimizer start
#   tables/numerical_checks.csv              dense-GLS cross-check of the engine
#   data/analysis_rows.csv.gz                exact rows entering every model
#   data/event_models.rds                    fitted models and designs
#                                            (input to scripts 09 and 10)
#   logs/R_session.txt                       R version and loaded packages
# =============================================================================

options(stringsAsFactors = FALSE)


# -----------------------------------------------------------------------------
# 0. Paths and shared functions
# -----------------------------------------------------------------------------

# Load 07_model_functions.R from this script's folder (Rscript), or from the
# working directory when the script is sourced interactively from the
# repository folder or from 5_did_cits/.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  getwd()
}
functions_file <- file.path(c(script_dir, file.path(script_dir, "5_did_cits")), "07_model_functions.R")
functions_file <- functions_file[file.exists(functions_file)][1]
if (is.na(functions_file)) stop("07_model_functions.R not found. Run this script from the repository folder.")
source(functions_file)   # paths, mixed_fit(), fit_checked(), effect_summary(), ...

repo_dir <- find_repository(dirname(functions_file))
data_dir <- repository_path(env_or("INDOOR_AQ_DATA", "data"), repo_dir)
args <- commandArgs(TRUE)
output_dir <- prepare_output(if (length(args) >= 2) argument_path(args[2]) else env_or("INDOOR_AQ_OUTPUT", "results"),
                             repo_dir, data_dir)
# Daily means: from step 03 in this run when available, otherwise the
# supplied data folder (see DATA_REQUIREMENTS.md).
input_dir <- if (length(args) >= 1) {
  argument_path(args[1])
} else if (file.exists(file.path(output_dir, "data/daily_home.csv"))) {
  file.path(output_dir, "data")
} else {
  data_dir
}
if (!file.exists(file.path(input_dir, "daily_home.csv"))) {
  stop("daily_home.csv not found. Run 03_prepare_cohort.ipynb first or supply it in data/; see DATA_REQUIREMENTS.md.")
}
for (file_name in c("events.csv", "home_id_crosswalk_23.csv")) {
  if (!file.exists(file.path(data_dir, file_name))) stop(paste(file_name, "not found in", data_dir, "; see DATA_REQUIREMENTS.md."))
}
cat("08_fit_did_cits_models.R\n  input :", input_dir, "\n  output:", output_dir, "\n")


# -----------------------------------------------------------------------------
# 1. Analysis settings
# -----------------------------------------------------------------------------

OUTCOMES <- c("temp", "rh", "pm25", "co2")
MODELS <- c("DiD", "CITS")
EFFECT_TERM <- "post:treated"           # retrofit effect in both models

PRE_DAYS <- 90                          # baseline window before works start
POST_DAYS <- 30                         # follow-up window after works end
MIN_PRE_DAYS <- 63                      # valid baseline days required per home
MIN_POST_DAYS <- 21                     # valid follow-up days required per home
MIN_COMPARISONS <- 5                    # comparison homes required per analysis
PM25_BASELINE_LIMIT <- 100              # ug/m3, primary PM2.5 analysis only
HOMES_WITH_OTHER_WORKS <- c("Home 4", "Home 17")  # never used as comparisons

HEATING_MONTHS <- c(11, 12, 1, 2, 3)    # November to March
TIME_UNIT_DAYS <- 30                    # time and posttime are in 30-day units
POSTTIME_CENTER <- 15.5                 # CITS effect at the follow-up midpoint
DAYNUM_ORIGIN <- as.Date("2024-01-01")  # day counter used for AR(1) gaps


# -----------------------------------------------------------------------------
# 2. Read the inputs
# -----------------------------------------------------------------------------

# Copy the three input files into the output run so that it is self-contained.
# Event dates and home identifiers come from the private data folder.
for (file_name in c("daily_home.csv", "events.csv", "home_id_crosswalk_23.csv")) {
  from <- if (file_name == "daily_home.csv") file.path(input_dir, file_name) else file.path(data_dir, file_name)
  to <- file.path(output_dir, "data", file_name)
  if (normalizePath(from, winslash = "/") != normalizePath(to, winslash = "/", mustWork = FALSE)) {
    if (!file.copy(from, to, overwrite = TRUE)) stop(paste("Could not copy input:", from))
  }
}

# Daily means (one row per home x outcome x local date; y = NA when the day
# had fewer than 18 valid hours).
daily <- read.csv(file.path(output_dir, "data/daily_home.csv"))
daily$date <- as.Date(daily$date)

# Retrofit events: treated home, first and last day of works.
events <- read.csv(file.path(output_dir, "data/events.csv"))
events$start <- as.Date(events$start)
events$end <- as.Date(events$end)

homes <- unique(daily$home)
treated_homes <- unique(events$home)
donors <- setdiff(homes, c(treated_homes, HOMES_WITH_OTHER_WORKS))  # comparison pool
stopifnot(length(homes) == 23, length(donors) == 17, nrow(events) == 7)


# -----------------------------------------------------------------------------
# 3. Functions for one event x outcome
#    `ev` is always one row of events.csv (event, home, start, end, ...).
# -----------------------------------------------------------------------------

bind_rows <- function(x) do.call(rbind, x)


#' All analysis dates: 90 baseline days and 30 follow-up days (works excluded).
analysis_dates <- function(ev) {
  c(seq(ev$start - PRE_DAYS, ev$start - 1, by = "day"),
    seq(ev$end + 1, ev$end + POST_DAYS, by = "day"))
}


#' Complete home x date grid for one outcome, with post and treated flags.
#' Dates without a valid daily mean keep y = NA.
event_outcome_grid <- function(ev, outcome) {
  grid <- merge(expand.grid(home = homes, date = analysis_dates(ev)),
                daily[daily$variable == outcome, c("home", "date", "y")],
                all.x = TRUE)
  grid$post <- as.integer(grid$date > ev$end)
  grid$treated <- as.integer(grid$home == ev$home)
  grid
}


#' Decide, for every home, whether it enters the primary analysis and the
#' expanded (SI) comparison pool, and record the reason.
screen_homes <- function(grid, ev, outcome) {
  bind_rows(lapply(homes, function(h) {
    q <- grid[grid$home == h, ]
    pre <- q$y[q$post == 0 & is.finite(q$y)]
    post <- q$y[q$post == 1 & is.finite(q$y)]

    candidate <- h == ev$home || h %in% donors
    coverage <- length(pre) >= MIN_PRE_DAYS && length(post) >= MIN_POST_DAYS
    baseline_ok <- outcome != "pm25" || (length(pre) > 0 && mean(pre) <= PM25_BASELINE_LIMIT)

    reason <- if (!candidate) {
      if (h %in% HOMES_WITH_OTHER_WORKS) "other_recorded_works" else "other_retrofitted_home"
    } else if (!coverage) {
      "insufficient_coverage"
    } else if (!baseline_ok) {
      "PM25_baseline_mean_above_100"
    } else {
      "included"
    }

    data.frame(
      event = ev$event, outcome = outcome, home = h,
      treated = h == ev$home, candidate = candidate,
      pre_n = length(pre), post_n = length(post),
      pre_mean = if (length(pre)) mean(pre) else NA,
      pre_median = if (length(pre)) median(pre) else NA,
      coverage_pass = coverage, baseline_pass = baseline_ok,
      primary_included = candidate && coverage && baseline_ok,
      original_pool_included = candidate && coverage,
      primary_reason = reason
    )
  }))
}


#' Observed before (pre) and after (post) means for every included home.
#' These feed Table 3 (treated home) and Table S5 (comparison homes).
observed_period_means <- function(z, selected, ev, outcome, analysis) {
  rows <- list()
  for (h in selected) {
    for (p in 0:1) {
      ix <- which(z$home == h & z$post == p)
      rows[[length(rows) + 1]] <- data.frame(
        event = ev$event, outcome = outcome, analysis = analysis, home = h,
        treated = h == ev$home, period = if (p) "post" else "pre",
        n = length(ix), mean = mean(z$y[ix]), median = median(z$y[ix])
      )
    }
  }
  bind_rows(rows)
}


#' Add the model variables (season, day counter, time, posttime, factors),
#' sort by home and day, and check the rows.
add_model_variables <- function(z, ev, analysis) {
  z$y_raw <- z$y                               # untouched copy, checked below
  z$event <- ev$event
  z$analysis <- analysis
  z$home <- factor(z$home)
  z$occasion <- factor(ev$event)
  z$block <- interaction(z$home, z$occasion, drop = TRUE)
  z$season <- factor(ifelse(as.integer(format(z$date, "%m")) %in% HEATING_MONTHS,
                            "heating", "non_heating"))
  z$daynum <- as.numeric(z$date - DAYNUM_ORIGIN)
  z$time <- as.numeric(z$date - ev$start) / TIME_UNIT_DAYS
  z$posttime <- z$post * (as.numeric(z$date - ev$end) - POSTTIME_CENTER) / TIME_UNIT_DAYS
  z <- z[order(z$home, z$daynum), ]
  rownames(z) <- NULL

  # One row per home-day, no construction days, response not transformed.
  stopifnot(!anyDuplicated(z[, c("home", "date")]),
            !any(z$date >= ev$start & z$date <= ev$end),
            identical(z$y, z$y_raw))
  z
}


#' Right-hand side of the DiD or CITS formula (season only if both occur).
model_rhs <- function(model, z) {
  rhs <- if (model == "DiD") "post*treated" else "post*treated+time*treated+posttime*treated"
  if (nlevels(z$season) > 1) rhs <- paste(rhs, "+season")
  rhs
}


#' Build the design for one model and fit it from three starting values.
fit_event_model <- function(z, model) {
  rhs <- model_rhs(model, z)
  X <- model.matrix(as.formula(paste("~", rhs)), z)
  stopifnot(qr(X)$rank == ncol(X))             # design must be full rank
  design <- list(z = z, X = X, Z = NULL, key = EFFECT_TERM, formula = paste("y ~", rhs))
  list(model = fit_checked(design), design = design)
}


#' One row of event_results.csv.
result_row <- function(ev, outcome, analysis, model, fit, selected) {
  m <- fit$model
  z <- fit$design$z
  effect <- effect_summary(m, EFFECT_TERM)
  data.frame(
    event = ev$event, home = ev$home, start = ev$start, end = ev$end,
    outcome = outcome, analysis = analysis, model = model,
    estimate = effect[["estimate"]], model_se = effect[["model_se"]],
    lower = effect[["lower"]], upper = effect[["upper"]], df = effect[["df"]],
    interval = "approximate_model_conditional_95pct",
    n_treated = 1, n_controls = m$G - 1, n_rows = m$n,
    n_unique_home_dates = nrow(unique(z[, c("home", "date")])),
    pre_n = sum(z$treated == 1 & z$post == 0),
    post_n = sum(z$treated == 1 & z$post == 1),
    winsorized = FALSE,
    PM25_baseline_limit = if (outcome == "pm25" && analysis == "primary") PM25_BASELINE_LIMIT else NA,
    phi = m$phi, var_home = m$var_home, var_resid = m$var_resid,
    objective = m$objective, convergence = m$convergence,
    formula = fit$design$formula,
    control_homes = paste(sort(setdiff(selected, ev$home)), collapse = ";")
  )
}


#' Independent check of the fast engine: rebuild the full covariance matrix
#' V explicitly, redo generalized least squares, and require the same
#' coefficients and covariance (to 1e-5).
check_against_dense_gls <- function(fit, ev, outcome, model) {
  m <- fit$model
  z <- fit$design$z
  X <- fit$design$X

  day_gap <- abs(outer(z$daynum, z$daynum, "-"))
  same_home <- outer(as.character(z$home), as.character(z$home), "==")
  V <- (m$var_resid * m$phi^day_gap + m$var_home) * same_home
  Vi <- chol2inv(chol(V))
  A <- solve(crossprod(X, Vi %*% X))
  beta_gls <- A %*% crossprod(X, Vi %*% z$y)

  beta_error <- max(abs(as.vector(beta_gls) - m$beta))
  vcov_error <- max(abs(A - m$vcov))
  stopifnot(beta_error < 1e-5, vcov_error < 1e-5)

  data.frame(event = ev$event, outcome = outcome, model = model,
             independent_GLS_beta_error = beta_error,
             independent_GLS_covariance_error = vcov_error,
             passed = TRUE)
}


# -----------------------------------------------------------------------------
# 4. Fit every event x outcome x analysis x model
# -----------------------------------------------------------------------------

results <- list()
eligibility <- list()
observed <- list()
analysis_rows <- list()
optimizer_starts <- list()
models <- list()
numerical_checks <- list()

for (i in seq_len(nrow(events))) {
  ev <- events[i, ]

  for (outcome in OUTCOMES) {
    grid <- event_outcome_grid(ev, outcome)
    selection <- screen_homes(grid, ev, outcome)
    eligibility[[length(eligibility) + 1]] <- selection

    # PM2.5 is also fitted with the expanded comparison pool (no baseline limit).
    analyses <- if (outcome == "pm25") c("primary", "SI_original_comparisons") else "primary"

    for (analysis in analyses) {
      included <- if (analysis == "primary") selection$primary_included else selection$original_pool_included
      selected <- selection$home[included]
      stopifnot(ev$home %in% selected, sum(selected %in% donors) >= MIN_COMPARISONS)

      # Rows with a valid daily mean for the selected homes.
      z <- grid[grid$home %in% selected & is.finite(grid$y), ]
      observed[[length(observed) + 1]] <- observed_period_means(z, selected, ev, outcome, analysis)
      z <- add_model_variables(z, ev, analysis)

      rows_used <- z
      rows_used$outcome <- outcome
      analysis_rows[[length(analysis_rows) + 1]] <- rows_used

      for (model in MODELS) {
        fit <- fit_event_model(z, model)
        m <- fit$model

        results[[length(results) + 1]] <- result_row(ev, outcome, analysis, model, fit, selected)

        audit <- m$start_audit
        audit$event <- ev$event
        audit$outcome <- outcome
        audit$model <- model
        audit$analysis <- analysis
        optimizer_starts[[length(optimizer_starts) + 1]] <- audit

        models[[paste(ev$event, outcome, analysis, model, sep = "_")]] <- fit

        # One representative pair of models is cross-checked by dense GLS.
        if (ev$event == "H16_Bundle" && outcome == "pm25" && analysis == "primary") {
          numerical_checks[[length(numerical_checks) + 1]] <- check_against_dense_gls(fit, ev, outcome, model)
        }

        cat(ev$event, outcome, analysis, model,
            sprintf("estimate=%.6f controls=%d", effect_summary(m)[["estimate"]], m$G - 1), "\n")
        flush.console()
      }
    }
  }

  # Save progress after every event (useful for long cluster jobs).
  write.csv(bind_rows(results), file.path(output_dir, "tables/event_results.csv"), row.names = FALSE)
}


# -----------------------------------------------------------------------------
# 5. Write the outputs
# -----------------------------------------------------------------------------

all_results <- bind_rows(results)
is_primary <- all_results$analysis == "primary"

write.csv(all_results[is_primary, ],
          file.path(output_dir, "tables/primary_event_results.csv"), row.names = FALSE)
write.csv(all_results[!is_primary, ],
          file.path(output_dir, "tables/SI_PM25_original_comparisons.csv"), row.names = FALSE)
write.csv(bind_rows(eligibility),
          file.path(output_dir, "tables/home_event_eligibility.csv"), row.names = FALSE, na = "")
write.csv(bind_rows(observed),
          file.path(output_dir, "tables/home_event_observed_means.csv"), row.names = FALSE)
write.csv(bind_rows(optimizer_starts),
          file.path(output_dir, "tables/optimizer_checks.csv"), row.names = FALSE)
write.csv(bind_rows(numerical_checks),
          file.path(output_dir, "tables/numerical_checks.csv"), row.names = FALSE)

connection <- gzfile(file.path(output_dir, "data/analysis_rows.csv.gz"), "wt")
write.csv(bind_rows(analysis_rows), connection, row.names = FALSE)
close(connection)

saveRDS(models, file.path(output_dir, "data/event_models.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs/R_session.txt"))

stopifnot(length(results) == 70, all(all_results$convergence == 0))
cat("COMPLETE: 56 primary and 14 SI event models; all unwinsorized; no pooled effect\n")
