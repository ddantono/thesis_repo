# methods/msvar/msvar.R
#
# Called by MATLAB via:
#   Rscript msvar.R <in_file> <cfg_file> <out_file>
#
# Dependencies (must be in the same directory):
#   msVAR_function.r   (Dallakyan et al. 2022, mmc2)
#   tsGLASSO.r         (Dallakyan et al. 2022, mmc3)
#
# Requires R packages: LongMemoryTS
#
# Output format (read by read_r_output_file.m):
#   Section "B"    : [K*p x K] coefficient matrix, row-major, lag-major order
#   Section "meta" : key,value pairs
#     Mandatory : lambda, selected, mse_per_eq, sparsity, n_nonzero
#     Optional  : p_selected, fdr_q, stage1_info, n_pairs_sig, oos_msfe

suppressPackageStartupMessages({
  if (!requireNamespace("LongMemoryTS", quietly = TRUE)) {
    stop("R package 'LongMemoryTS' is required. Install with: install.packages('LongMemoryTS')")
  }
  library(LongMemoryTS)
})

# ── Command-line arguments ────────────────────────────────────────────────────
# args: <in_file> <cfg_file> <script_dir> <out_file>
# script_dir passed explicitly by MATLAB — avoids sys.frame() issues
args       <- commandArgs(trailingOnly = TRUE)
in_file    <- args[1]
cfg_file   <- args[2]
script_dir <- args[3]   # passed as-is from MATLAB; already a valid path
out_file   <- args[4]

# ── Source companion scripts ──────────────────────────────────────────────────
setwd(script_dir)
source(file.path(script_dir, "tsGLASSO.r"))
source(file.path(script_dir, "msVAR_function.r"))

# ── Read input ────────────────────────────────────────────────────────────────
Y   <- as.matrix(read.csv(in_file))
cfg <- read.csv(cfg_file, stringsAsFactors = FALSE)
cfg <- setNames(as.list(cfg$value), cfg$param)

K                  <- as.integer(cfg$K)
p_seq              <- as.integer(scan(text = cfg$p_seq_str, quiet = TRUE))
lambda_cfg         <- as.numeric(cfg$lambda)
half_window_length <- as.integer(cfg$half_window_length)
rho                <- as.numeric(cfg$rho)
alpha_admm         <- as.numeric(cfg$alpha)
thresh             <- as.numeric(cfg$thresh)
rho_flex           <- as.logical(as.integer(cfg$rho_flex))
ADMM_iter          <- as.integer(cfg$ADMM_iter)
stage1_info        <- cfg$stage1_info          # 'bic' or 'aic'
standardize        <- as.logical(as.integer(cfg$standardize))
lambda_gam         <- as.numeric(cfg$lambda_gam)
lambda_trim_max    <- as.numeric(cfg$lambda_trim_max)
lambda_trim_min    <- as.numeric(cfg$lambda_trim_min)
lambda_criteria    <- cfg$lambda_criteria      # 'BIC','AIC','eBIC','CV'
fdr_q              <- as.numeric(cfg$fdr_q)
ite_max            <- as.integer(cfg$ite_max)
reltol             <- as.numeric(cfg$reltol)

stopifnot(ncol(Y) == K)
stopifnot(length(p_seq) >= 1, all(p_seq >= 0))

T_obs <- nrow(Y)
p_ub  <- max(p_seq)

# ── Lambda selection ──────────────────────────────────────────────────────────
# sentinel -1 means auto-select via BIC over a lambda grid
if (is.na(lambda_cfg) || lambda_cfg < 0) {
  lambda <- tryCatch({
    selectlambda(X                = Y,
                 gam              = lambda_gam,
                 thresh           = thresh,
                 halfWindowLength = half_window_length,
                 rho              = rho,
                 Max_Iter         = ADMM_iter,
                 diag             = FALSE,
                 verbose          = FALSE,
                 rho.flex         = rho_flex,
                 alpha            = alpha_admm,
                 trim_max         = lambda_trim_max,
                 trim_min         = lambda_trim_min,
                 freq             = NULL,
                 smooth           = FALSE,
                 criteria         = lambda_criteria)$lambda
  }, error = function(e) {
    # Fallback: analytic upper bound
    2 * sqrt(log(p_ub) / T_obs)
  })
} else {
  lambda <- lambda_cfg
}

# Scalar guard (eBIC/all criteria return a vector)
if (length(lambda) > 1) {
  lambda <- lambda[1]
}

# ── Run msVAR ─────────────────────────────────────────────────────────────────
# Strategy:
#   1. Try with TSGlasso spectral screening (standard msVAR)
#   2. If TSGlasso finds no significant pairs (subscript OOB), retry with
#      allPairs = all K*(K-1)/2 pairs — equivalent to sVAR Stage 1 with
#      full connectivity, letting Stage 2 FDR do the selection.
#   3. If both fail, stop with informative message.

i_run_msvar <- function(all_pairs_override = NULL) {
  tryCatch(
    msVAR(dta                  = Y,
          p.seq                = p_seq,
          lambda               = lambda,
          halfWindowLength     = half_window_length,
          stage1.showStatus    = FALSE,
          stage2.showStatus    = FALSE,
          standardize          = standardize,
          stage1.info          = stage1_info,
          iteMax               = ite_max,
          reltol               = reltol,
          rho.flex             = rho_flex,
          rho                  = rho,
          alpha                = alpha_admm,
          fdr.q                = fdr_q,
          ADMM_ITER            = ADMM_iter,
          allPairs             = all_pairs_override),
    error = function(e) {
      message("msVAR() ERROR: ", conditionMessage(e))
      NULL
    }
  )
}

# Attempt 1: standard spectral screening
msvar_result <- i_run_msvar(all_pairs_override = NULL)

# Attempt 2: if failed (TSGlasso found no pairs), use all pairs
if (is.null(msvar_result)) {
  message("msVAR attempt 1 failed — retrying with allPairs (full connectivity).")
  all_pairs_mat <- do.call(rbind, unlist(
    lapply(1:(K-1), function(i) lapply((i+1):K, function(j) c(i,j))),
    recursive = FALSE))
  msvar_result <- i_run_msvar(all_pairs_override = all_pairs_mat)
}

if (is.null(msvar_result)) {
  stop("msVAR() failed on both attempts — check data (T, K) and R packages.")
}

# ── Extract Stage-2 output ───────────────────────────────────────────────────
# Primary output: Stage 2 (FDR-refined); fall back to Stage 1 if Stage 2 is NULL
stage2 <- msvar_result$stage2
stage1 <- msvar_result$stage1

if (is.null(stage2)) {
  # p == 0 path: model collapsed to white noise; use zero matrix
  A_est   <- matrix(0, K, K * p_ub)
  p_sel   <- 0L
  warning("msVAR Stage 2 is NULL (p=0 model selected). Returning zero B matrix.")
} else {
  A_est <- as.matrix(stage2$estA)   # [K x K*p_sel], row-indexed by equation
  p_sel <- as.integer(stage2$order)
}

# ── Convert A from [K x K*p] (row=eq, col=var_lag) to B [K*p x K] ───────────
# R convention:  A_est[i, (k-1)*K + j]  = coefficient of y_j(t-k) in eq i
# MATLAB/framework convention: B[(k-1)*K + j, i]  (lag-major, same meaning)
# => B = t(A_est)
if (is.null(stage2) || p_sel == 0) {
  B <- matrix(0, K * p_ub, K)
} else {
  # A_est may have fewer lag columns than K*p_ub if p_sel < p_ub; pad with zeros
  if (ncol(A_est) < K * p_ub) {
    A_est <- cbind(A_est, matrix(0, K, K * p_ub - ncol(A_est)))
  }
  B <- t(A_est)   # [K*p_ub x K]
}

# ── Compute fitted values, residuals, diagnostics ────────────────────────────
build_lagged <- function(Y, p, K) {
  T_obs <- nrow(Y)
  Y0    <- Y[(p + 1):T_obs, , drop = FALSE]
  X     <- matrix(0, T_obs - p, K * p)
  for (lag in seq_len(p)) {
    cols      <- ((lag - 1) * K + 1):(lag * K)
    X[, cols] <- Y[(p + 1 - lag):(T_obs - lag), , drop = FALSE]
  }
  list(Y0 = Y0, X = X)
}

if (p_sel == 0L) {
  # White-noise model: Y_t = intercept + u_t
  intercept <- if (!is.null(stage1) && !is.null(stage1$estIntercept)) {
    as.numeric(stage1$estIntercept)
  } else {
    colMeans(Y)
  }
  
  residuals  <- sweep(Y, 2, intercept, "-")
  mse_per_eq <- colMeans(residuals^2)
  
} else {
  lm_obj <- build_lagged(Y, p_sel, K)
  Y0     <- lm_obj$Y0
  X      <- lm_obj$X
  
  B_fit <- B[1:(K * p_sel), , drop = FALSE]
  
  intercept <- if (!is.null(stage2$estIntercept)) {
    as.numeric(stage2$estIntercept)
  } else {
    rep(0, K)
  }
  
  fitted    <- sweep(X %*% B_fit, 2, intercept, "+")
  residuals <- Y0 - fitted
  
  mse_per_eq <- colMeans(residuals^2)
}

sparsity  <- mean(B == 0)
n_nonzero <- sum(B != 0)
n_pairs   <- if (!is.null(msvar_result$spectral)) {
  tryCatch(nrow(msvar_result$spectral$gpre.max), error = function(e) NA)
} else NA

# ── Write output ──────────────────────────────────────────────────────────────
con <- file(out_file, "w")

# Section B: [K*p_ub x K] coefficient matrix
writeLines("B", con)
for (row in seq_len(nrow(B))) {
  writeLines(paste(B[row, ], collapse = ","), con)
}

# Section meta: key,value pairs
writeLines("meta", con)
# Write intercept vector so MATLAB can assess centering quality
if (!is.null(stage2) && !is.null(stage2$estIntercept)) {
  intc <- as.numeric(stage2$estIntercept)
  writeLines(paste0("intercept_max,",  max(abs(intc))),            con)
  writeLines(paste0("intercept_mean,", mean(abs(intc))),           con)
  writeLines(paste0("intercept_centered,",
                    as.integer(max(abs(intc)) < 1e-6)),             con)
} else if (!is.null(stage1) && !is.null(stage1$estIntercept)) {
  intc <- as.numeric(stage1$estIntercept)
  writeLines(paste0("intercept_max,",  max(abs(intc))),            con)
  writeLines(paste0("intercept_mean,", mean(abs(intc))),           con)
  writeLines(paste0("intercept_centered,",
                    as.integer(max(abs(intc)) < 1e-6)),             con)
}
writeLines(paste0("lambda,",      lambda),                         con)
writeLines(paste0("selected,", paste(as.integer(t(B != 0)), collapse = ",")), con)
writeLines(paste0("mse_per_eq,",  paste(mse_per_eq, collapse = ",")), con)
writeLines(paste0("sparsity,",    sparsity),                       con)
writeLines(paste0("n_nonzero,",   n_nonzero),                      con)
writeLines(paste0("p_selected,",  p_sel),                          con)
writeLines(paste0("p_ub,",        p_ub),                           con)
writeLines(paste0("fdr_q,",       fdr_q),                          con)
writeLines(paste0("stage1_info,", stage1_info),                    con)
writeLines(paste0("oos_msfe,",    NaN),                            con)  # no rolling CV
if (!is.na(n_pairs)) {
  writeLines(paste0("n_pairs_sig,", n_pairs),                      con)
}

close(con)