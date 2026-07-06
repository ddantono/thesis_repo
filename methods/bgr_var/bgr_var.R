# methods/bgr_var/bgr_var.R
#
# BGR-VAR via BigVAR — called by run_bgr_var.m via system call.
#
# USAGE (from MATLAB call_r_script):
#   Rscript bgr_var.R <in_file> <cfg_file> <out_file>
#
# INPUT FILES
#   in_file  : CSV with header y1..yK — the Y matrix [T x K]
#   cfg_file : CSV with two columns param,value — scalar/string config
#
# OUTPUT FILE
#   out_file : sections B and meta (read by read_r_output_file.m)
#
# OUTPUT FORMAT (must match read_r_output_file.m conventions)
#   Section "B"    : K*p rows x K cols, comma-separated doubles
#   Section "meta" : key,value pairs
#     Mandatory : lambda, mse, sparsity, n_nonzero, selected
#     Optional  : oos_msfe, struct
#
# Reference: Banbura, Giannone & Reichlin (2010)
#            Nicholson, Matteson & Bien (2017), BigVAR R package

suppressPackageStartupMessages({
  library(BigVAR)
})

# --------------------------------------------------------------------------- #
#  0. Parse command-line arguments                                             #
# --------------------------------------------------------------------------- #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: bgr_var.R <in_file> <cfg_file> <out_file>")
}

in_file  <- args[1]
cfg_file <- args[2]
out_file <- args[3]

# --------------------------------------------------------------------------- #
#  1. Read Y                                                                   #
# --------------------------------------------------------------------------- #
Y <- as.matrix(read.csv(in_file, header = TRUE))
T_obs <- nrow(Y)
K     <- ncol(Y)

# --------------------------------------------------------------------------- #
#  2. Read config                                                               #
# --------------------------------------------------------------------------- #
cfg_raw <- read.csv(cfg_file, header = TRUE, stringsAsFactors = FALSE)
cfg     <- setNames(as.list(cfg_raw$value), cfg_raw$param)

i_cfg_num <- function(key, default) {
  if (!is.null(cfg[[key]])) as.numeric(cfg[[key]]) else default
}
i_cfg_bool <- function(key, default) {
  if (!is.null(cfg[[key]])) as.logical(as.integer(cfg[[key]])) else default
}

p         <- as.integer(i_cfg_num("p",        1))
nlambda   <- as.integer(i_cfg_num("nlambda",  50))
n_folds   <- as.integer(i_cfg_num("n_folds",  10))
T1        <- as.integer(i_cfg_num("T1",       floor(T_obs / 3)))
T2        <- as.integer(i_cfg_num("T2",       floor(2 * T_obs / 3)))
Minnesota <- i_cfg_bool("Minnesota", FALSE)
intercept <- i_cfg_bool("intercept", TRUE)

# --------------------------------------------------------------------------- #
#  3. Build and fit BigVAR model (BGR structure)                               #
# --------------------------------------------------------------------------- #
# BGR = Bayesian Ridge Regression (Banbura et al. 2010).
# gran = c(nlambda, n_folds): grid depth and rolling CV evaluations.
# IC = FALSE: skip AIC/BIC benchmarks for speed.
# RVAR = FALSE: BGR is a ridge estimator — no sparse support to refit on.

model <- tryCatch({
  constructModel(
    Y,
    p       = p,
    struct  = "BGR",
    gran    = c(nlambda, n_folds),
    verbose = FALSE,
    IC      = FALSE,
    T1      = T1,
    T2      = T2,
    model.controls = list(
      MN        = Minnesota,
      intercept = intercept,
      RVAR      = FALSE
    )
  )
}, error = function(e) {
  stop(paste("constructModel failed:", conditionMessage(e)))
})

results <- tryCatch({
  cv.BigVAR(model)
}, error = function(e) {
  stop(paste("cv.BigVAR failed:", conditionMessage(e)))
})

# --------------------------------------------------------------------------- #
#  4. Extract coefficient matrix B                                             #
# --------------------------------------------------------------------------- #
# betaPred: k x (k*p + 1) — last column is intercept when intercept=TRUE.
# Framework convention: B is [K*p x K] (lag-major, no intercept column).
# BigVAR stores: rows = equations (k), cols = lag terms + intercept.
# We transpose and drop the intercept column.

beta_full <- results@betaPred   # [K x K*p+1] or [K x K*p]

# Drop intercept column (first column) if present
if (ncol(beta_full) == K * p + 1) {
  beta_no_intercept <- beta_full[, 2:(K * p + 1), drop = FALSE]
} else {
  beta_no_intercept <- beta_full
}

# Transpose to [K*p x K] — framework convention
B <- t(beta_no_intercept)   # [K*p x K]

# --------------------------------------------------------------------------- #
#  5. Compute in-sample metrics from betaPred                                  #
# --------------------------------------------------------------------------- #
# fitted and resids slots are computed by BigVAR from betaPred.
# Use them directly if available, otherwise recompute.

if (length(results@fitted) > 0 && length(results@resids) > 0) {
  fitted_vals <- results@fitted   # [T_eff x K] or similar
  resids_vals <- results@resids
} else {
  # Manual recomputation using betaPred and last-window Zvals
  # Fallback: compute from B and Y directly (lag-major)
  T_eff <- T_obs - p
  Y0    <- Y[(p + 1):T_obs, , drop = FALSE]
  X     <- matrix(0, nrow = T_eff, ncol = K * p)
  for (lag in 1:p) {
    cols         <- ((lag - 1) * K + 1):(lag * K)
    X[, cols]    <- Y[(p + 1 - lag):(T_obs - lag), , drop = FALSE]
  }
  fitted_vals <- X %*% B
  resids_vals <- Y0 - fitted_vals
}

# MSE per equation (in-sample, from residuals)
mse_per_eq <- colMeans(resids_vals^2)

# --------------------------------------------------------------------------- #
#  6. Compute summary statistics                                               #
# --------------------------------------------------------------------------- #
n_total   <- length(B)
# BGR is ridge — no exact zeros expected, but apply a small threshold
# for sparsity reporting (informational only, not for support recovery).
tol       <- 1e-6
n_nonzero <- sum(abs(B) > tol)

# selected: number of "active" coefficients above threshold
selected <- as.integer(abs(B) > tol)
sparsity  <- mean(selected == 0)

# OOS MSFE from BigVAR rolling CV — OOSMSFE is a vector (one per equation); take mean
oos_msfe_raw <- results@OOSMSFE
oos_msfe     <- if (length(oos_msfe_raw) > 0) mean(as.numeric(oos_msfe_raw), na.rm = TRUE) else NA_real_

# Optimal lambda — scalar
opt_lambda_raw <- results@OptimalLambda
opt_lambda     <- if (length(opt_lambda_raw) > 0) as.numeric(opt_lambda_raw)[1] else NA_real_

# --------------------------------------------------------------------------- #
#  7. Write output file                                                         #
# --------------------------------------------------------------------------- #
fcon <- file(out_file, "w")

# Section B: K*p rows x K cols
writeLines("B", fcon)
for (i in 1:nrow(B)) {
  writeLines(paste(sprintf("%.15g", B[i, ]), collapse = ","), fcon)
}

# Section meta: key,value pairs
writeLines("meta", fcon)
writeLines(paste("lambda",    sprintf("%.15g", opt_lambda), sep = ","), fcon)
writeLines(paste("mse",       paste(sprintf("%.15g", mse_per_eq), collapse = ","), sep = ","), fcon)
writeLines(paste("sparsity",  sprintf("%.15g", sparsity),  sep = ","), fcon)
writeLines(paste("n_nonzero", sprintf("%d",    n_nonzero), sep = ","), fcon)
writeLines(paste("selected",  paste(selected, collapse = ","), sep = ","), fcon)
writeLines(paste("oos_msfe",  sprintf("%.15g", if (!is.na(oos_msfe)) oos_msfe else 0.0), sep = ","), fcon)
writeLines(paste("struct",    "BGR",                                                          sep = ","), fcon)
writeLines(paste("nlambda",   sprintf("%d", nlambda),                                         sep = ","), fcon)
writeLines(paste("n_folds",   sprintf("%d", n_folds),                                         sep = ","), fcon)

close(fcon)

cat(sprintf("[bgr_var.R] Done. K=%d p=%d T=%d lambda=%.4g oos_msfe=%.6f\n",
            K, p, T_obs,
            if (!is.na(opt_lambda)) opt_lambda else 0.0,
            if (!is.na(oos_msfe))  oos_msfe  else 0.0))