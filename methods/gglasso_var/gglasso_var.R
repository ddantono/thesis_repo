# methods/gglasso_var/gglasso_var.R
#
# Group LASSO VAR estimation via the gglasso R package.
# Called by run_gglasso_var.m through call_r_script.m.
#
# ARGS (commandArgs):
#   args[1] : input CSV  (Y matrix, headers y1..yK)
#   args[2] : config CSV (param,value pairs)
#   args[3] : output CSV (written by this script)
#
# OUTPUT FORMAT (sections separated by bare section-name lines):
#   B
#     K*p rows x K columns — coefficient matrix
#   meta
#     key,value pairs (see REQUIRED keys below)
#
# REQUIRED meta keys:
#   lambda        : CV-selected lambda per equation (K values, comma-sep)
#   selected      : number of non-zero coefficients
#   mse_per_eq    : in-sample MSE per equation (K values)
#   sparsity      : fraction of zero coefficients in B
#   n_nonzero     : number of non-zero coefficients
#   oos_msfe      : NaN (gglasso CV is in-sample; no rolling OOS)
#   cv_criterion  : lambda selection criterion used
#   n_groups      : number of groups per equation
#   group_by      : grouping scheme used

suppressPackageStartupMessages({
  library(gglasso)
})

# --------------------------------------------------------------------------
# 0. Parse arguments
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || any(nchar(args[1:3]) == 0)) {
  stop("[gglasso_var.R] Usage: Rscript gglasso_var.R <in> <cfg> <out>")
}
in_file  <- trimws(args[1])
cfg_file <- trimws(args[2])
out_file <- trimws(args[3])

# --------------------------------------------------------------------------
# 1. Read data
# --------------------------------------------------------------------------
Y <- as.matrix(read.csv(in_file, header = TRUE, check.names = FALSE))
Y <- apply(Y, 2, as.numeric)   # force numeric
T_obs <- nrow(Y)
K     <- ncol(Y)

# --------------------------------------------------------------------------
# 2. Read config
# --------------------------------------------------------------------------
cfg_raw <- read.csv(cfg_file, header = TRUE, stringsAsFactors = FALSE)
cfg     <- setNames(as.list(cfg_raw$value), cfg_raw$param)

p            <- as.integer(cfg$p); if(is.na(p)) stop("p not specified in cfg")
nlambda      <- as.integer(cfg$nlambda); if(is.na(nlambda)) nlambda <- 100
loss         <- as.character(cfg$loss); if(is.na(loss)) loss <- "ls"
cv_criterion <- as.character(cfg$cv_criterion); if(is.na(cv_criterion)) cv_criterion <- "lambda.min"
n_folds      <- as.integer(cfg$n_folds); if(is.na(n_folds)) n_folds <- 5
eps_tol      <- as.numeric(cfg$eps); if(is.na(eps_tol)) eps_tol <- 1e-6
maxit        <- as.integer(cfg$maxit); if(is.na(maxit)) maxit <- 1000
use_intercept <- as.logical(as.integer(cfg$intercept)); if(is.na(use_intercept)) use_intercept <- TRUE
group_by     <- as.character(cfg$group_by); if(is.na(group_by)) group_by <- "lag"
lambda_user  <- as.numeric(cfg$lambda); if(is.na(lambda_user)) lambda_user <- -1

# --------------------------------------------------------------------------
# 3. Build lag-major regressor matrix X and response matrix Y0
#    X columns: [y1(t-1)..yK(t-1), y1(t-2)..yK(t-2), ..., y1(t-p)..yK(t-p)]
# --------------------------------------------------------------------------
T_eff <- T_obs - p
if (T_eff <= 0) stop("Number of observations T_obs <= p; cannot construct lagged matrix")
Y0    <- Y[(p + 1):T_obs, , drop = FALSE]   # [T_eff x K]

X <- matrix(0, nrow = T_eff, ncol = K * p)
for (lag in 1:p) {
  cols       <- ((lag - 1) * K + 1):(lag * K)
  row_from   <- (p + 1 - lag)
  row_to     <- (T_obs - lag)
  X[, cols]  <- Y[row_from:row_to, , drop = FALSE]
}
n_pred <- ncol(X)   # K*p

# --------------------------------------------------------------------------
# 3b. Center Y0 and X when intercept=FALSE (remove column means)
#     gglasso with intercept=FALSE does NOT center internally.
#     Without centering, non-zero column means bias the coefficients.
# --------------------------------------------------------------------------
if (!use_intercept) {
  Y0_means <- colMeans(Y0)
  X_means  <- colMeans(X)
  Y0 <- sweep(Y0, 2, Y0_means, "-")
  X  <- sweep(X, 2, X_means, "-")
}

# --------------------------------------------------------------------------
# 4. Build group index vector
#    "lag"      : all K predictors at the same lag form one group
#                 group vector: rep(1:p, each=K)  → K*p entries, p groups/eq
#    "variable" : all p lags of the same predictor form one group
#                 group vector: rep(1:K, times=p) → K*p entries, K groups/eq
# --------------------------------------------------------------------------
if (group_by == "variable") {
  group_vec <- rep(1:K, times = p)
  n_groups  <- K
} else if (group_by == "lag") {
  group_vec <- rep(1:p, each = K)
  n_groups  <- p
} else {
  warning(sprintf("Unknown group_by='%s', defaulting to 'lag'", group_by))
  group_vec <- rep(1:p, each = K)
  n_groups  <- p
}

# --------------------------------------------------------------------------
# 5. Fit group lasso equation-by-equation
# --------------------------------------------------------------------------
B          <- matrix(0, nrow = n_pred, ncol = K)
lambda_sel <- numeric(K)
mse_per_eq <- numeric(K)
intercept_vals <- numeric(K)   # store intercept coefficients per equation

for (k in 1:K) {
  y_k <- Y0[, k]
  
  tryCatch({
    if (lambda_user < 0) {
      # Auto-select via CV
      cv_fit <- cv.gglasso(
        x       = X,
        y       = y_k,
        group   = group_vec,
        loss    = loss,
        pred.loss = "L2",
        nlambda = nlambda,
        nfolds  = n_folds,
        eps     = eps_tol,
        maxit   = maxit,
        intercept = use_intercept
      )
      lam_use <- if (cv_criterion == "lambda.1se") {
        cv_fit$lambda.1se
      } else {
        cv_fit$lambda.min
      }
      fit <- cv_fit$gglasso.fit
    } else {
      # User-supplied lambda — fit single model (no CV)
      fit    <- gglasso(
        x       = X,
        y       = y_k,
        group   = group_vec,
        loss    = loss,
        lambda  = lambda_user,
        eps     = eps_tol,
        maxit   = maxit,
        intercept = use_intercept
      )
      lam_use <- lambda_user
    }
    
    # Extract coefficients at selected lambda
    # coef() returns (n_pred+1) x 1 including intercept in row 1
    coef_full  <- as.numeric(coef(fit, s = lam_use))
    intercept_vals[k] <- coef_full[1] 
    # Drop intercept (row 1); keep predictors
    B[, k]     <- coef_full[-1]
    lambda_sel[k] <- lam_use
    
  }, error = function(e) {
    warning(sprintf("[gglasso_var.R] Equation %d failed: %s", k, e$message))
    # B[, k] stays zero; lambda_sel[k] stays 0
  })
  
  # In-sample MSE from MATLAB-side residuals uses the MATLAB B,
  # but we compute here for the meta section (R-side, consistent).
  fitted_k    <- X %*% B[, k]
  resid_k     <- y_k - fitted_k
  mse_per_eq[k] <- mean(resid_k^2)
}

# --------------------------------------------------------------------------
# 6. Compute global diagnostics
# --------------------------------------------------------------------------
sparsity <- sum(B == 0, na.rm = TRUE) / length(B)
n_nonzero <- sum(B != 0, na.rm = TRUE)
fitted_vals <- sweep(X %*% B, 2, intercept_vals, "+")
residuals   <- Y0 - fitted_vals
mse         <- colMeans(residuals^2)

# --------------------------------------------------------------------------
# 7. Write output
# --------------------------------------------------------------------------
fout <- file(out_file, "w", encoding = "UTF-8")

writeLines("B", fout)
for (i in 1:nrow(B)) {
  writeLines(paste(B[i, ], collapse = ","), fout)
}

writeLines("meta", fout)
writeLines(paste0("lambda,",       paste(lambda_sel, collapse = ",")), fout)
writeLines(paste0("selected,",     paste(as.integer(B != 0), collapse = ",")), fout)
writeLines(paste0("mse,",          paste(round(mse, 15), collapse = ",")), fout)
writeLines(paste0("sparsity,",     round(sparsity, 15)), fout)
writeLines(paste0("n_nonzero,",    sum(B != 0)), fout)
writeLines(paste0("oos_msfe,",     "NaN"), fout)
writeLines(paste0("cv_criterion,", cv_criterion), fout)
writeLines(paste0("n_groups,",     n_groups), fout)
writeLines(paste0("group_by,",     group_by), fout)
writeLines(paste0("intercept,",    paste(intercept_vals, collapse = ",")), fout)

close(fout)