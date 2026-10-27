# =============================================================================
# 10_run_placebo_tests.R  (step 10)
#
# Treated-label placebo (permutation-style) tests.
#
# Manuscript : Table S6 (56 primary models) and Results 3.3 (Home 14 CWI CITS
#              PM2.5 model with the expanded comparison pool)
# Usage      : Rscript 5_did_cits/10_run_placebo_tests.R [MODEL_RUN] [OUTPUT_DIR]
#                MODEL_RUN  folder containing data/event_models.rds from
#                           08_fit_did_cits_models.R. Default: OUTPUT_DIR
#                OUTPUT_DIR results folder. Default: results/
#                           (INDOOR_AQ_OUTPUT)
# Runtime    : one to two minutes (~600 refits x 3 optimizer starts)
#
# Method
# ------
# For each saved model, keep exactly the same rows (homes, dates, outcome
# values) and formula, but pretend in turn that each participating home was
# the treated one, rebuild the design matrix and refit. The participating
# homes are never reselected. The observed assignment is included, so with G
# homes there are G refits.
#
#   rank = number of finite refits whose |estimate| >= |observed estimate|
#   p    = rank / number of finite refits   (two-sided empirical tail)
#
# Refits that fail are recorded (convergence = 999, estimate = NA) and are
# excluded from the denominator. No placebo is run on the log scale.
#
# Outputs (inside OUTPUT_DIR)
# ---------------------------
#   tables/placebo_assignments.csv          every refit (one row per home)
#   tables/placebo_summary.csv              rank and p for all 57 models
#   tables/placebo_ranks_all_56_primary.csv the 56 primary models (Table S6)
#   logs/placebo_R_session.txt              R version and loaded packages
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
cat("10_run_placebo_tests.R\n  input :", input_dir, "\n  output:", output_dir, "\n")


# -----------------------------------------------------------------------------
# 1. Models to test: the 56 primary models plus one expanded-pool model
# -----------------------------------------------------------------------------

EXPANDED_POOL_KEY <- "H14_CWI_pm25_SI_original_comparisons_CITS"

saved <- readRDS(file.path(input_dir, "data/event_models.rds"))
keys <- c(grep("_primary_", names(saved), value = TRUE), EXPANDED_POOL_KEY)
stopifnot(length(keys) == 57, all(keys %in% names(saved)))


# -----------------------------------------------------------------------------
# 2. Functions
# -----------------------------------------------------------------------------

#' Copy a stored design, label home `h` as the treated home, and rebuild X.
#' Rows, dates, homes and outcome values are checked to be unchanged.
relabel_treated_home <- function(design, h) {
  z <- design$z
  placebo <- design
  placebo$z$treated <- as.integer(as.character(z$home) == h)
  placebo$X <- model.matrix(as.formula(design$formula), placebo$z)

  stopifnot(identical(placebo$z$y, z$y),
            identical(placebo$z$date, z$date),
            identical(placebo$z$home, z$home),
            qr(placebo$X)$rank == ncol(placebo$X))
  placebo
}


#' Refit one model once per participating home (G refits).
#' @return One row per assigned home with its estimate and convergence code.
placebo_refits <- function(key, design) {
  z <- design$z
  observed_home <- unique(as.character(z$home[z$treated == 1]))
  stopifnot(length(observed_home) == 1)
  homes <- sort(unique(as.character(z$home)))

  rows <- list()
  for (h in homes) {
    placebo <- relabel_treated_home(design, h)
    # The observed assignment must reproduce the original design exactly.
    if (h == observed_home) stopifnot(identical(placebo$X, design$X))

    m <- tryCatch(fit_checked(placebo), error = function(e) e)
    failed <- inherits(m, "error")
    rows[[h]] <- data.frame(
      key = key,
      assigned_treated = h,
      observed_assignment = h == observed_home,
      estimate = if (failed) NA_real_ else unname(m$beta[design$key]),
      convergence = if (failed) 999 else m$convergence,
      error = if (failed) conditionMessage(m) else ""
    )
  }
  do.call(rbind, rows)
}


#' Empirical two-sided rank and tail probability of the observed estimate.
placebo_summary <- function(key, design, refits) {
  observed <- refits$estimate[refits$observed_assignment]
  stopifnot(length(observed) == 1, is.finite(observed))

  finite <- is.finite(refits$estimate)
  rank <- sum(abs(refits$estimate[finite]) >= abs(observed))
  info <- model_key_info(key, design)

  data.frame(
    event = info$event, outcome = info$outcome,
    analysis = info$analysis, model = info$model,
    estimate = observed,
    G = nrow(refits),
    finite_refits = sum(finite),
    rank = rank,
    p = rank / sum(finite),
    failed_refits = sum(!finite)
  )
}


# -----------------------------------------------------------------------------
# 3. Run all placebo distributions
# -----------------------------------------------------------------------------

assignments <- list()
summaries <- list()

for (key in keys) {
  design <- saved[[key]]$design
  refits <- placebo_refits(key, design)
  summary_row <- placebo_summary(key, design, refits)

  assignments[[key]] <- refits
  summaries[[key]] <- summary_row

  # Save progress after every model.
  write.csv(do.call(rbind, assignments),
            file.path(output_dir, "tables/placebo_assignments.csv"), row.names = FALSE)
  write.csv(do.call(rbind, summaries),
            file.path(output_dir, "tables/placebo_summary.csv"), row.names = FALSE)

  cat(summary_row$event, summary_row$outcome, summary_row$model, summary_row$analysis,
      "rank", summary_row$rank, "of", summary_row$finite_refits, "\n")
  flush.console()
}


# -----------------------------------------------------------------------------
# 4. Write the Table S6 input
# -----------------------------------------------------------------------------

all_summaries <- do.call(rbind, summaries)
write.csv(all_summaries[all_summaries$analysis == "primary", ],
          file.path(output_dir, "tables/placebo_ranks_all_56_primary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs/placebo_R_session.txt"))

cat("Complete: 56 primary placebo distributions and the reported expanded-pool Home 14 CWI CITS check.\n")
