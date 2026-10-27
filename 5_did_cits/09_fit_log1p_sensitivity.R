# =============================================================================
# 09_fit_log1p_sensitivity.R  (step 09)
#
# PM2.5 log-scale sensitivity analysis: refit the 14 primary PM2.5 models
# (7 events x DiD/CITS) after replacing the response y with log(1 + y).
#
# Manuscript : Table S3 (log1p columns) and Supplementary Methods S3
#              (residual excess kurtosis)
# Usage      : Rscript 5_did_cits/09_fit_log1p_sensitivity.R [MODEL_RUN] [OUTPUT_DIR]
#                MODEL_RUN  folder containing data/event_models.rds from
#                           08_fit_did_cits_models.R. Default: OUTPUT_DIR
#                OUTPUT_DIR results folder. Default: results/
#                           (INDOOR_AQ_OUTPUT)
# Runtime    : a few seconds
#
# Only the response changes. Each refit reuses the stored primary design
# exactly: same homes, dates, covariates, design matrix and eligibility
# (eligibility stays based on the original-scale baseline mean). This is
# checked before every fit. Log-scale coefficients are differences on the
# log(1 + y) scale, not concentration differences. No placebo analysis is
# run on the log scale.
#
# Outputs (inside OUTPUT_DIR)
# ---------------------------
#   tables/paired_PM25_results.csv   original vs log1p estimates, intervals
#                                    and residual excess kurtosis
#   tables/log_design_checks.csv     confirms identical rows and design
#   tables/log_optimizer_starts.csv  outcome of every optimizer start
#   data/matched_log1p_models.rds    fitted log1p models and designs
#   logs/log_R_session.txt           R version and loaded packages
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
# Saved primary models from step 08 (by default in the same results folder).
input_dir <- if (length(args) >= 1) argument_path(args[1]) else output_dir
if (!file.exists(file.path(input_dir, "data/event_models.rds"))) {
  stop("data/event_models.rds not found. Run 08_fit_did_cits_models.R first.")
}
cat("09_fit_log1p_sensitivity.R\n  input :", input_dir, "\n  output:", output_dir, "\n")


# -----------------------------------------------------------------------------
# 1. Select the 14 primary PM2.5 models
# -----------------------------------------------------------------------------

saved <- readRDS(file.path(input_dir, "data/event_models.rds"))
keys <- grep("_pm25_primary_", names(saved), value = TRUE)
stopifnot(length(keys) == 14)


# -----------------------------------------------------------------------------
# 2. Refit each model on the log(1 + y) scale
# -----------------------------------------------------------------------------

#' Copy a stored design and replace only the response with log(1 + y).
log1p_design <- function(design) {
  # The stored response must be the untransformed, non-negative concentration.
  stopifnot(is.null(design$Z),
            all(is.finite(design$z$y)),
            all(design$z$y >= 0),
            identical(design$z$y, design$z$y_raw))

  log_design <- design
  log_design$z$y <- log1p(design$z$y)

  # Everything except y must be unchanged.
  other_columns <- setdiff(names(design$z), "y")
  stopifnot(identical(log_design$X, design$X),
            identical(log_design$Z, design$Z),
            identical(log_design$z[, other_columns], design$z[, other_columns]))
  log_design
}

results <- list()
models <- list()
design_checks <- list()
optimizer_starts <- list()

for (key in keys) {
  original <- saved[[key]]
  log_design <- log1p_design(original$design)
  log_model <- fit_checked(log_design)

  raw <- effect_summary(original$model)
  transformed <- effect_summary(log_model)
  event <- as.character(unique(original$design$z$event))
  model <- sub(".*_primary_", "", key)

  results[[key]] <- data.frame(
    event = event, model = model,
    raw_estimate = raw[["estimate"]],
    raw_lower = raw[["lower"]],
    raw_upper = raw[["upper"]],
    log1p_estimate = transformed[["estimate"]],
    log1p_lower = transformed[["lower"]],
    log1p_upper = transformed[["upper"]],
    raw_residual_excess_kurtosis = residual_kurtosis(original$model$residuals),
    log1p_residual_excess_kurtosis = residual_kurtosis(log_model$residuals),
    same_sign = sign(raw[["estimate"]]) == sign(transformed[["estimate"]]),
    log1p_convergence = log_model$convergence,
    n_controls = log_model$G - 1,
    n_rows = log_model$n
  )
  models[[key]] <- list(model = log_model, design = log_design)
  design_checks[[key]] <- data.frame(event = event, model = model,
                                     rows_dates_covariates_identical = TRUE,
                                     fixed_design_identical = TRUE)

  audit <- log_model$start_audit
  audit$event <- event
  audit$model <- model
  optimizer_starts[[key]] <- audit

  cat(event, model, "log1p estimate", transformed[["estimate"]], "\n")
  flush.console()

  # Save progress after every model.
  write.csv(do.call(rbind, results), file.path(output_dir, "tables/paired_PM25_results.csv"), row.names = FALSE)
}


# -----------------------------------------------------------------------------
# 3. Write the outputs
# -----------------------------------------------------------------------------

saveRDS(models, file.path(output_dir, "data/matched_log1p_models.rds"))
write.csv(do.call(rbind, design_checks),
          file.path(output_dir, "tables/log_design_checks.csv"), row.names = FALSE)
write.csv(do.call(rbind, optimizer_starts),
          file.path(output_dir, "tables/log_optimizer_starts.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs/log_R_session.txt"))

cat("Complete: 14 matched log1p fits. No placebo analysis is performed on the log scale.\n")
