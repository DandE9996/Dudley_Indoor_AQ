# =============================================================================
# 07_model_functions.R  (step 07)
#
# Functions shared by the R scripts 08_fit_did_cits_models.R,
# 09_fit_log1p_sensitivity.R and 10_run_placebo_tests.R. Base R only.
# This file only defines functions; it reads no data and fits no model.
#
#   1. Repository paths  find_repository(), repository_path(), env_or(),
#                        argument_path(), prepare_output()
#   2. Model engine      mixed_fit(): REML linear mixed model with a home
#                        random intercept and AR(1) errors
#   3. Model helpers     fit_checked(), effect_summary(), residual_kurtosis(),
#                        model_key_info()
# =============================================================================


# =============================================================================
# 1. Repository paths
#
# The R scripts use the same folders as the notebooks: supplied data in data/
# (INDOOR_AQ_DATA) and generated results in results/ (INDOOR_AQ_OUTPUT).
# Relative paths are resolved against the repository folder, not the shell.
# =============================================================================

#' Repository folder: INDOOR_AQ_CODE, or the first folder at or above `start`
#' that contains 00_shared_functions.ipynb and the 5_did_cits/ folder.
find_repository <- function(start) {
  configured <- Sys.getenv("INDOOR_AQ_CODE")
  if (nzchar(configured)) return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  dir <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "00_shared_functions.ipynb")) &&
        dir.exists(file.path(dir, "5_did_cits"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) stop("Repository folder not found; set INDOOR_AQ_CODE.")
    dir <- parent
  }
}


#' Value of an environment variable, or `default` when it is unset or empty.
env_or <- function(name, default) {
  value <- Sys.getenv(name)
  if (nzchar(value)) value else default
}


#' Absolute path; relative paths are taken from the repository folder.
repository_path <- function(value, repo_dir) {
  value <- gsub("\\\\", "/", path.expand(value))
  absolute <- startsWith(value, "/") || grepl("^[A-Za-z]:/", value)
  if (!absolute) value <- file.path(repo_dir, value)
  normalizePath(value, winslash = "/", mustWork = FALSE)
}


#' Absolute path for a command-line argument; relative paths are taken from
#' the current working directory, as usual for command-line tools.
argument_path <- function(value) {
  normalizePath(path.expand(value), winslash = "/", mustWork = FALSE)
}


#' Create the results folder (data/, tables/, logs/) and refuse the
#' repository folder itself and the supplied data folder.
prepare_output <- function(value, repo_dir, data_dir) {
  out <- repository_path(value, repo_dir)
  key <- function(x) if (.Platform$OS.type == "windows") tolower(x) else x
  protected <- data_dir
  inside <- vapply(protected, function(p) {
    key(out) == key(p) || startsWith(key(out), paste0(key(p), "/"))
  }, logical(1))
  if (key(out) == key(repo_dir) || any(inside)) {
    stop("Choose a results folder separate from the supplied data folder.")
  }
  for (sub_dir in c("data", "tables", "logs")) {
    dir.create(file.path(out, sub_dir), recursive = TRUE, showWarnings = FALSE)
  }
  out
}


# =============================================================================
# 2. Model engine
#
# Linear mixed model used for every DiD and CITS fit in the manuscript:
#
#     y = X beta + b_home + e
#
#   * X beta : fixed effects (DiD or CITS design, optional heating season)
#   * b_home : home random intercept,  b_home ~ N(0, lam * sig)
#   * e      : AR(1) errors within each home, using the ACTUAL number of
#              calendar days between consecutive observations, so the
#              correlation between two days that are d days apart is phi^d
#              (missing days therefore weaken the correlation correctly).
#
# The model is fitted by exact Gaussian restricted maximum likelihood (REML).
# The residual variance `sig` and the fixed effects `beta` are profiled out
# analytically, so the optimizer only searches over
#
#     par[1] = log(lam)   variance ratio  var(home intercept) / var(residual)
#     par[2] = phi        AR(1) correlation for a one-day gap, in [0, 0.995]
#     par[3] = log(tau)   smoothing-variance ratio (only when a spline basis Z
#                         is supplied; the manuscript models never use one)
#
# How the likelihood is evaluated (function `calc` below)
# -------------------------------------------------------
# 1. Prewhitening. Within each home, row t is replaced by
#        (row_t - rho_t * row_{t-1}) / sqrt(1 - rho_t^2),  rho_t = phi^gap_t,
#    which turns irregularly spaced AR(1) errors into independent errors.
#    The first row of each home is left unchanged.
# 2. Random intercept. After prewhitening, the column of ones becomes the
#    vector `v`. The home intercept adds lam * v v' to the error covariance,
#    which is inverted in closed form with the Sherman-Morrison formula
#    (scalar `adj`) and whose log-determinant comes from the matrix
#    determinant lemma. No n x n matrix is ever built.
# 3. Each home contributes a small cross-product matrix of [X, Z, y].
#    Summing them gives X'V^-1 X, X'V^-1 y and y'V^-1 y, from which the GLS
#    estimate, the residual sum of squares and the REML objective follow.
#
# Implementation notes
# --------------------
# * y is standardized (mean ym, SD ys) before fitting for numerical stability;
#   all returned coefficients, variances and residuals are back on the
#   original scale. X is not standardized; column 1 must be the intercept.
# * The returned `vcov` is the model-based covariance used for every interval
#   in the manuscript. `cluster_vcov` (a home-clustered sandwich estimator) is
#   returned for completeness only and is not used in the paper.
# * `residuals` are conditional residuals: observed y minus the fixed-effect
#   fit minus the predicted home intercept.
#
# The code of mixed_fit() is unchanged from the analysis reported in the
# article; only comments and layout were edited.
# =============================================================================

#' Fit the REML mixed model described above.
#'
#' @param ds    Model specification, a list with elements
#'                z       data frame with columns y, home, occasion, daynum,
#'                        date and block (one row per home-day)
#'                X       fixed-effect design matrix (column 1 = intercept)
#'                Z       optional penalized spline basis, or NULL
#'                formula text of the model formula (stored in the result)
#' @param ar    TRUE to estimate the AR(1) correlation; FALSE fixes phi = 0.
#' @param start Starting values c(log(lam), phi[, log(tau)]); NULL uses
#'              c(0, 0.5[, 0]).
#' @param maxit Maximum number of L-BFGS-B iterations.
#' @return A list with the fixed effects (`beta`), their covariance (`vcov`),
#'         variance components, residuals, the optimizer convergence code
#'         (`convergence`, 0 = converged) and the REML objective (`objective`,
#'         -2 x restricted log-likelihood up to a constant; smaller is better).
mixed_fit <- function(ds, ar = TRUE, start = NULL, maxit = 100) {

  # ---------------------------------------------------------------------------
  # Sort the rows by home, occasion and day, and standardize the response.
  #   n = rows, p = fixed effects, k = spline coefficients (0 here),
  #   M = [X, Z, y], G = number of homes (clusters)
  # ---------------------------------------------------------------------------
  z <- ds$z
  ix <- order(z$home, z$occasion, z$daynum)
  z <- z[ix, ]
  X <- ds$X[ix, , drop = FALSE]
  Z <- if (!is.null(ds$Z)) ds$Z[ix, , drop = FALSE] else matrix(numeric(), nrow(X), 0)

  n <- nrow(X)
  p <- ncol(X)
  k <- ncol(Z)
  ym <- mean(z$y)
  ys <- sd(z$y)
  y <- (z$y - ym) / ys
  M <- cbind(X, Z, y)
  indices <- split(seq_len(n), z$home, drop = TRUE)
  G <- length(indices)

  # For each home: its rows of M, the day gaps between consecutive rows, and
  # whether a new occasion (event) starts, which resets the AR(1) chain.
  blocks <- lapply(indices, function(i) {
    a <- z[i, ]
    delta <- diff(a$daynum)
    new <- as.character(a$occasion[-1]) != as.character(a$occasion[-nrow(a)])
    list(M = M[i, , drop = FALSE], delta = delta, new = new)
  })

  # ---------------------------------------------------------------------------
  # REML objective for one parameter vector.
  #   detail = FALSE : return the objective only (used by the optimizer)
  #   detail = TRUE  : return the full fitted model
  # A non-positive-definite configuration returns 1e60 (or NULL).
  # ---------------------------------------------------------------------------
  calc <- function(par, detail = FALSE) {
    lam <- exp(par[1])                    # home-intercept variance ratio
    phi <- if (ar) par[2] else 0          # AR(1) correlation per day
    tau <- if (k) exp(par[3]) else 0      # spline variance ratio
    S <- matrix(0, ncol(M), ncol(M))      # running sum of M' V^-1 M
    ld <- 0                               # running log-determinant of V
    cache <- list()

    for (j in seq_along(blocks)) {
      b <- blocks[[j]]
      a <- b$M
      nr <- nrow(a)

      # Step 1: AR(1) prewhitening with gap-dependent correlation phi^gap.
      rho <- phi^pmax(b$delta, 1)
      rho[b$new] <- 0
      den <- sqrt(1 - rho^2)
      W <- a
      if (nr > 1) W[-1, ] <- (a[-1, , drop = FALSE] - rho * a[-nr, , drop = FALSE]) / den

      # Step 2: random intercept via Sherman-Morrison.
      #   v   = prewhitened column of ones
      #   adj = lam / (1 + lam * v'v)
      v <- c(1, (1 - rho) / den)
      u <- as.vector(crossprod(v, W))
      adj <- lam / (1 + lam * sum(v^2))

      # Step 3: this home's contribution to M' V^-1 M and to log|V|.
      Sj <- crossprod(W) - adj * tcrossprod(u)
      S <- S + Sj
      ld <- ld + sum(log1p(-rho^2)) + log1p(lam * sum(v^2))
      if (detail) cache[[j]] <- list(S = Sj, W = W, v = v, adj = adj, u = u)
    }

    # Keep the fixed-effect columns and y (drop the spline columns).
    keep <- c(seq_len(p), ncol(M))
    A0 <- S[keep, keep, drop = FALSE]
    B <- NULL
    Kinv <- NULL

    # Optional penalized spline: integrate out its coefficients (Schur
    # complement). Not used by the manuscript models (k = 0).
    if (k) {
      zi <- p + seq_len(k)
      K <- S[zi, zi, drop = FALSE] + diag(1 / tau, k)
      chK <- tryCatch(chol(K), error = function(e) NULL)
      if (is.null(chK)) return(if (detail) NULL else 1e60)
      Kinv <- chol2inv(chK)
      B <- S[keep, zi, drop = FALSE]
      A0 <- A0 - B %*% Kinv %*% t(B)
      ld <- ld + k * log(tau) + 2 * sum(log(diag(chK)))
    }

    # Generalized least squares for the fixed effects (standardized scale):
    #   A  = X'V^-1 X,   bb = A^-1 X'V^-1 y,   rss = y'V^-1 y - bb'X'V^-1 y
    A <- A0[seq_len(p), seq_len(p), drop = FALSE]
    ch <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(ch)) return(if (detail) NULL else 1e60)

    Ai <- chol2inv(ch)
    bb <- as.vector(Ai %*% A0[seq_len(p), p + 1])
    rss <- A0[p + 1, p + 1] - sum(bb * A0[seq_len(p), p + 1])
    if (!is.finite(rss) || rss <= 0) return(if (detail) NULL else 1e60)

    # Profiled REML objective (-2 x restricted log-likelihood + constant):
    #   (n - p) log(sig) + log|V| + log|X'V^-1 X|
    sig <- rss / (n - p)
    obj <- (n - p) * log(sig) + ld + 2 * sum(log(diag(ch)))
    if (!detail) return(obj)

    # ---- Everything below runs only once, at the optimum (detail = TRUE) ----

    # Predicted spline coefficients (empty when k = 0).
    ufit <- if (k) as.vector(Kinv %*% (S[zi, ncol(M)] - S[zi, seq_len(p), drop = FALSE] %*% bb)) else numeric()

    # Per-home score contributions X_j' V_j^-1 e_j for the sandwich estimator.
    score <- matrix(0, G, p)
    for (j in seq_len(G)) {
      a <- cache[[j]]
      er <- a$W[, ncol(M)] - as.vector(a$W[, seq_len(p), drop = FALSE] %*% bb)
      if (k) er <- er - as.vector(a$W[, zi, drop = FALSE] %*% ufit)
      wx <- a$W[, seq_len(p), drop = FALSE]
      score[j, ] <- as.vector(crossprod(wx, er) - a$adj * a$u[seq_len(p)] * sum(a$v * er))
    }

    # Home-clustered sandwich covariance with small-sample correction,
    # rescaled to the original units (reported for completeness only).
    cv <- Ai %*% crossprod(score) %*% Ai * G / (G - 1) * (n - 1) / (n - p) * ys^2

    # Fixed effects on the original scale (the mean of y returns to the intercept).
    beta <- bb * ys
    beta[1] <- beta[1] + ym
    names(beta) <- colnames(X)

    # Conditional residuals: subtract the fixed-effect fit and each home's
    # predicted (BLUP) intercept `bhat`.
    fitfixed <- as.vector(X %*% beta + if (k) Z %*% (ufit * ys) else 0)
    res <- z$y - fitfixed
    for (j in seq_len(G)) {
      a <- cache[[j]]
      r <- a$W[, ncol(M)] - as.vector(a$W[, seq_len(p), drop = FALSE] %*% bb)
      if (k) r <- r - as.vector(a$W[, zi, drop = FALSE] %*% ufit)
      bhat <- a$adj * sum(a$v * r) * ys
      res[indices[[j]]] <- res[indices[[j]]] - bhat
    }

    list(
      beta         = beta,                # fixed effects
      vcov         = Ai * sig * ys^2,     # model-based covariance of beta (used in the paper)
      cluster_vcov = cv,                  # sandwich covariance (not used in the paper)
      smooth_coeff = ufit * ys,
      var_home     = lam * sig * ys^2,    # home random-intercept variance
      var_resid    = sig * ys^2,          # marginal AR(1) residual variance
      var_smooth   = tau * sig * ys^2,
      phi          = phi,                 # AR(1) correlation for a one-day gap
      n            = n,                   # rows (home-days)
      G            = G,                   # homes (clusters)
      p            = p,                   # fixed effects
      objective    = obj,                 # REML objective at the optimum
      residuals    = res,                 # conditional residuals
      fitted       = z$y - res,
      fixed_fitted = fitfixed,
      dates        = z$date,
      home         = z$home,
      block        = z$block,
      X            = X,
      smooth_Z     = Z,
      formula      = ds$formula
    )
  }

  # ---------------------------------------------------------------------------
  # Optimize the REML objective with box-constrained L-BFGS-B.
  # ---------------------------------------------------------------------------
  if (is.null(start)) start <- c(0, .5, if (k) 0)
  if (!ar) start[2] <- 0

  if (ar) {
    # Search log(lam) in [-16, 16], phi in [0, 0.995] (and log(tau) in [-16, 16]).
    op <- optim(start, calc, method = 'L-BFGS-B',
                lower = c(-16, 0, if (k) -16),
                upper = c(16, .995, if (k) 16),
                control = list(maxit = maxit, factr = 1e8))
  } else {
    # Without AR(1), keep phi fixed at 0 and optimize the remaining parameters.
    free <- c(1, if (k) 3)
    op <- optim(start[free],
                function(a) {
                  pp <- start
                  pp[free] <- a
                  pp[2] <- 0
                  calc(pp)
                },
                method = 'L-BFGS-B',
                lower = rep(-16, length(free)),
                upper = rep(16, length(free)),
                control = list(maxit = maxit, factr = 1e8))
    pp <- start
    pp[free] <- op$par
    pp[2] <- 0
    op$par <- pp
  }

  # Evaluate the full model at the optimum and record the optimizer status.
  ans <- calc(op$par, TRUE)
  if (is.null(ans)) stop('marginal_likelihood_failed')
  ans$convergence <- op$convergence
  ans$parameters <- op$par
  ans
}


# =============================================================================
# 3. Model helpers
# =============================================================================

# Starting values c(log(lam), phi) tried for every model: the engine default
# (log(lam) = 0, phi = 0.5), then a weak and a strong home/AR(1) structure.
OPTIMIZER_STARTS <- list(NULL, c(-3, .15), c(3, .8))


#' Fit a model from several starting values and keep the best solution.
#'
#' Each start is fitted with mixed_fit(). Among the starts that converged,
#' the one with the smallest REML objective is returned. A start that fails
#' with an error is recorded with convergence code 999.
#'
#' @param ds Model specification passed to mixed_fit().
#' @return The best fitted model, with an extra element `start_audit`: one row
#'         per start with its convergence code, REML objective and the
#'         `post:treated` estimate (lets reviewers see that starts agree).
fit_checked <- function(ds) {
  fits <- lapply(OPTIMIZER_STARTS, function(s) {
    tryCatch(mixed_fit(ds, start = s, maxit = 250), error = function(e) NULL)
  })

  converged <- vapply(fits, function(m) !is.null(m) && m$convergence == 0, logical(1))
  if (!any(converged)) stop("No converged optimization")

  candidates <- fits[converged]
  objectives <- vapply(candidates, `[[`, numeric(1), "objective")
  best <- candidates[[which.min(objectives)]]

  best$start_audit <- do.call(rbind, lapply(seq_along(fits), function(j) {
    m <- fits[[j]]
    if (is.null(m)) {
      return(data.frame(start = j, convergence = 999, objective = NA, estimate = NA))
    }
    data.frame(start = j, convergence = m$convergence, objective = m$objective,
               estimate = unname(m$beta["post:treated"]))
  }))
  best
}


#' Point estimate and approximate model-conditional 95% interval.
#'
#' The interval uses the model-based standard error and a t distribution with
#' G - 2 degrees of freedom (G = number of homes in the model), because each
#' event has a single treated home.
#'
#' @param m   Fitted model returned by mixed_fit() / fit_checked().
#' @param key Name of the coefficient; `post:treated` is the retrofit effect
#'            in both the DiD and the CITS design.
#' @return Named numeric vector: estimate, model_se, lower, upper, df.
effect_summary <- function(m, key = "post:treated") {
  j <- match(key, names(m$beta))
  estimate <- unname(m$beta[j])
  se <- sqrt(m$vcov[j, j])
  df <- m$G - 2
  c(estimate = estimate,
    model_se = se,
    lower    = estimate - qt(.975, df) * se,
    upper    = estimate + qt(.975, df) * se,
    df       = df)
}


#' Residual excess kurtosis (Supplementary Methods S3).
#'
#' Fourth central moment divided by the fourth power of the sample standard
#' deviation, minus 3. A normal distribution gives approximately 0.
residual_kurtosis <- function(x) {
  mean((x - mean(x))^4) / sd(x)^4 - 3
}


#' Describe one saved model from its key in event_models.rds.
#'
#' Keys have the form <event>_<outcome>_<analysis>_<model>, for example
#' "H14_CWI_pm25_SI_original_comparisons_CITS". Event and analysis are read
#' from the stored design rows because both can contain underscores.
#'
#' @param key    Name of the element in event_models.rds.
#' @param design The stored design (`$design`) of that element.
#' @return A list with event, outcome, analysis and model.
model_key_info <- function(key, design) {
  event <- as.character(unique(design$z$event))
  outcome <- sub(paste0("^", event, "_"), "", key)
  outcome <- sub("_.*", "", outcome)
  list(event = event,
       outcome = outcome,
       analysis = as.character(unique(design$z$analysis)),
       model = if (grepl("_CITS$", key)) "CITS" else "DiD")
}
