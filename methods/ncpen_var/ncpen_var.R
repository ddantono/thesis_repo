# methods/ncpen_var/ncpen_var.R
#
# Called by MATLAB via:
#   Rscript ncpen_var.R <in_file> <cfg_file> <out_file>
#
# Requires: BigVAR

suppressPackageStartupMessages(library(BigVAR))

args <- commandArgs(trailingOnly = TRUE)
in_file  <- args[1]
cfg_file <- args[2]
out_file <- args[3]

# --- Read input ---------------------------------------------------------
Y   <- as.matrix(read.csv(in_file))
cfg <- read.csv(cfg_file, stringsAsFactors = FALSE)
cfg <- setNames(as.list(cfg$value), cfg$param)

p         <- as.integer(cfg$p)
K_cfg     <- as.integer(cfg$K)
penalty   <- cfg$penalty          # expected: "SCAD" or "MCP"
n_folds   <- as.integer(cfg$n_folds)
nlambda   <- as.integer(cfg$nlambda)
Minnesota <- as.logical(as.integer(cfg$Minnesota))
verbose   <- as.logical(as.integer(cfg$verbose))

K <- ncol(Y)
if (!is.na(K_cfg) && K_cfg != K) {
  stop(sprintf("Config K=%d does not match ncol(Y)=%d", K_cfg, K))
}

if (!(penalty %in% c("SCAD", "MCP"))) {
  stop("penalty must be 'SCAD' or 'MCP'")
}

# BigVAR uses gran = c(grid_depth, number_of_cv_evaluations)
gran <- c(nlambda, n_folds)

# --- Fit BigVAR model ---------------------------------------------------
mod <- constructModel(
  Y,
  p       = p,
  struct  = penalty,
  gran    = gran,
  h       = 1,
  cv      = "Rolling",
  verbose = verbose,
  IC      = FALSE,
  model.controls = list(
    intercept = TRUE,
    MN        = Minnesota,
    gamma     = 3
  )
)

res <- cv.BigVAR(mod)

# --- Extract coefficient matrix -----------------------------------------
# betaPred is K x (K*p + 1) for plain VAR with intercept:
# col 1 = intercept, cols 2:(K*p+1) = lag coefficients
beta <- res@betaPred
has_intercept <- (ncol(beta) == K * p + 1)

if (has_intercept) {
  intercept <- as.numeric(beta[, 1])
  B_slopes  <- beta[, 2:(K * p + 1), drop = FALSE]   # K x (K*p)
} else {
  intercept <- rep(0, K)
  B_slopes  <- beta
}

# Framework convention: B is (K*p) x K
B <- t(B_slopes)

# --- Build lag-major X and compute fitted values consistently ----------
build_lagged <- function(Y, p, K) {
  T_obs <- nrow(Y)
  Y0    <- Y[(p + 1):T_obs, , drop = FALSE]
  X     <- matrix(0, nrow = T_obs - p, ncol = K * p)

  for (lag in 1:p) {
    cols      <- ((lag - 1) * K + 1):(lag * K)
    X[, cols] <- Y[(p + 1 - lag):(T_obs - lag), , drop = FALSE]
  }

  list(Y0 = Y0, X = X)
}

lm_obj <- build_lagged(Y, p, K)
Y0     <- lm_obj$Y0
X      <- lm_obj$X

fitted_vals <- sweep(X %*% B, 2, intercept, "+")
residuals   <- Y0 - fitted_vals
mse         <- colMeans(residuals^2)

# --- Diagnostics --------------------------------------------------------
lambda_opt_raw <- res@OptimalLambda
lambda_opt     <- if (length(lambda_opt_raw) > 0) as.numeric(lambda_opt_raw)[1] else NA_real_
sparse_count <- res@sparse_count
selected     <- as.integer(abs(B) > 1e-8)
n_nonzero    <- sum(selected)
sparsity     <- mean(selected == 0)
oos_msfe     <- mean(as.numeric(res@OOSMSFE), na.rm = TRUE)

lag_active <- matrix(0, nrow = p, ncol = K)
for (eq in 1:K) {
  for (lag in 1:p) {
    rows <- ((lag - 1) * K + 1):(lag * K)
    lag_active[lag, eq] <- as.integer(any(abs(B[rows, eq]) > 1e-8))
  }
}

# --- Write output -------------------------------------------------------
con <- file(out_file, "w")

writeLines("B", con)
for (row in 1:nrow(B)) {
  writeLines(paste(B[row, ], collapse = ","), con)
}

writeLines("meta", con)
writeLines(paste0("lambda,",       lambda_opt), con)
writeLines(paste0("selected,",     paste(selected, collapse = ",")), con)
writeLines(paste0("mse,",          paste(mse, collapse = ",")), con)
writeLines(paste0("sparsity,",     sparsity), con)
writeLines(paste0("n_nonzero,",    n_nonzero), con)
writeLines(paste0("cv_sparsity,",  1 - sparse_count), con)
writeLines(paste0("oos_msfe,",     oos_msfe), con)
writeLines(paste0("penalty,",      penalty), con)
writeLines(paste0("n_folds,",      n_folds), con)
writeLines(paste0("nlambda,",      nlambda), con)
writeLines(paste0("gran_depth,",   nlambda), con)
writeLines(paste0("cv_type,",      "Rolling"), con)
writeLines(paste0("intercept,",    paste(intercept, collapse = ",")), con)
writeLines(paste0("lag_active,",   paste(as.integer(lag_active), collapse = ",")), con)

close(con)