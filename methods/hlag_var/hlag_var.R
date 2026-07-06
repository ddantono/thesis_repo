# methods/hlag_var/hlag_var.R
#
# Called by MATLAB via:
#   Rscript hlag_var.R <in_file> <cfg_file> <out_file>
#
# Supports: struct = 'HLAGOO' or 'HLAGC'
# Requires: BigVAR

suppressPackageStartupMessages(library(BigVAR))

args     <- commandArgs(trailingOnly = TRUE)
in_file  <- args[1]
cfg_file <- args[2]
out_file <- args[3]

# --- Read input ---------------------------------------------------------
Y   <- as.matrix(read.csv(in_file))
cfg <- read.csv(cfg_file, stringsAsFactors = FALSE)
cfg <- setNames(as.list(cfg$value), cfg$param)

p             <- as.integer(cfg$p)
K             <- as.integer(cfg$K)
struct        <- cfg$struct
nlambda       <- as.integer(cfg$nlambda)
n_folds       <- as.integer(cfg$n_folds)
T1            <- as.integer(cfg$T1)
T2            <- as.integer(cfg$T2)
RVAR          <- as.logical(as.integer(cfg$RVAR))
Minnesota     <- as.logical(as.integer(cfg$Minnesota))
intercept     <- as.logical(as.integer(cfg$intercept))
rel_threshold <- if (!is.null(cfg$rel_threshold)) as.numeric(cfg$rel_threshold) else 0.0  # 0 = no thresholding

if (!(struct %in% c('HLAGOO', 'HLAGC'))) {
  stop(paste('Invalid struct:', struct))
}

T_obs <- nrow(Y)
T_eff <- T_obs - p

# --- Construct and cross-validate BigVAR model -------------------------
model <- constructModel(
  Y,
  p       = p,
  struct  = struct,
  gran    = c(nlambda, n_folds),
  T1      = T1,
  T2      = T2,
  verbose = FALSE,
  IC      = FALSE,
  model.controls = list(
    RVAR      = RVAR,
    MN        = Minnesota,
    intercept = intercept
  )
)

results    <- cv.BigVAR(model)
lambda_opt <- results@OptimalLambda
oos_msfe   <- mean(as.numeric(results@OOSMSFE), na.rm = TRUE)

# --- Extract betaPred ---------------------------------------------------
# betaPred: coefficient matrix from the last rolling window.
# # Standard BigVAR output for prediction purposes.
B_raw <- results@betaPred                      # [K x K*p+1]
B_raw <- B_raw[, 2:(K * p + 1), drop = FALSE] # drop intercept column
B     <- t(B_raw)                             # [K*p x K]

# --- Optional relative threshold ----------------------------------------
# rel_threshold = 0.0  → no post-processing (raw BigVAR output)
# rel_threshold = 0.05 → zero entries with |b| < 5% of max(|B|)
#
# Note: HLAG does not produce exact zeros — it produces very small values
# for non-selected terms due to the hierarchical group penalty structure.
# The appropriate threshold depends on the intended use:
#   - Forecasting (OOS MSFE): use rel_threshold = 0.0
#   - Sparsity visualization: use rel_threshold = 0.05
#   - Support recovery: HLAG is not designed for this task
if (rel_threshold > 0) {
  max_b <- max(abs(B))
  if (max_b > 0) {
    B[abs(B) < rel_threshold * max_b] <- 0
  }
}

# --- Fitted values and residuals ----------------------------------------
Y0 <- Y[(p + 1):T_obs, , drop = FALSE]
X  <- matrix(0, nrow = T_eff, ncol = K * p)
for (lag in 1:p) {
  cols      <- ((lag - 1) * K + 1):(lag * K)
  X[, cols] <- Y[(p + 1 - lag):(T_obs - lag), , drop = FALSE]
}

fitted    <- X %*% B
residuals <- Y0 - fitted
mse       <- colMeans(residuals^2)
selected  <- as.integer(B != 0)
n_nonzero <- sum(selected)
sparsity  <- mean(selected == 0)

# --- Write output -------------------------------------------------------
con <- file(out_file, 'w')

writeLines('B', con)
for (row in 1:nrow(B)) {
  writeLines(paste(B[row, ], collapse = ','), con)
}

writeLines('meta', con)
writeLines(paste0('lambda,',        lambda_opt),                      con)
writeLines(paste0('selected,',      paste(selected, collapse = ',')),  con)
writeLines(paste0('mse,',           paste(mse, collapse = ',')),      con)
writeLines(paste0('sparsity,',      sparsity),                        con)
writeLines(paste0('n_nonzero,',     n_nonzero),                       con)
writeLines(paste0('oos_msfe,',      oos_msfe),                        con)
writeLines(paste0('struct,',        struct),                          con)
writeLines(paste0('nlambda,',       nlambda),                         con)
writeLines(paste0('n_folds,',       n_folds),                         con)
writeLines(paste0('rel_threshold,', rel_threshold),                   con)

close(con)