# methods/adaptive_lasso_var/adaptive_lasso_var.R
#
# Called by MATLAB via:
#   Rscript adaptive_lasso_var.R <in_file> <cfg_file> <out_file>
#
# Requires: glmnet

suppressPackageStartupMessages(library(glmnet))

args     <- commandArgs(trailingOnly = TRUE)
in_file  <- args[1]
cfg_file <- args[2]
out_file <- args[3]

# --- Read input ---------------------------------------------------------
Y   <- as.matrix(read.csv(in_file))
cfg <- read.csv(cfg_file, stringsAsFactors = FALSE)
cfg <- setNames(as.list(cfg$value), cfg$param)

p            <- as.integer(cfg$p)
K            <- as.integer(cfg$K)
nu           <- as.numeric(cfg$nu)
weight_eps   <- as.numeric(cfg$weight_eps)
n_folds      <- as.integer(cfg$n_folds)
nlambda      <- as.integer(cfg$nlambda)
standardize  <- as.logical(as.integer(cfg$standardize))
alpha        <- as.numeric(cfg$alpha)
cv_criterion <- cfg$cv_criterion   # 'lambda.min' or 'lambda.1se'
pilot        <- cfg$pilot          # 'ols' or 'lasso'

T_obs <- nrow(Y)
T_eff <- T_obs - p

# --- Build lag-major regressor matrix X ---------------------------------
build_lagged <- function(Y, p, K) {
  T_obs <- nrow(Y)
  Y0    <- Y[(p+1):T_obs, , drop=FALSE]
  X     <- matrix(0, nrow=T_obs-p, ncol=K*p)
  for (lag in 1:p) {
    cols        <- ((lag-1)*K + 1):(lag*K)
    X[, cols]   <- Y[(p+1-lag):(T_obs-lag), , drop=FALSE]
  }
  list(Y0=Y0, X=X)
}

lm_obj  <- build_lagged(Y, p, K)
Y0      <- lm_obj$Y0
X       <- lm_obj$X

# --- Equation-by-equation Adaptive LASSO --------------------------------
B          <- matrix(0, nrow=K*p, ncol=K)
lambda_sel <- numeric(K)

pilot_cfg <- pilot    # FIX: save configured pilot BEFORE loop

for (k in 1:K) {
  y_k   <- Y0[, k]
  pilot <- pilot_cfg  # reset to configured value each iteration
  
  # Stage 1: pilot
  if (pilot == 'ols') {
    b_pilot <- tryCatch(
      as.numeric(coef(lm(y_k ~ X - 1))),
      error = function(e) NULL
    )
    if (is.null(b_pilot) || any(!is.finite(b_pilot))) {
      pilot <- 'lasso'   # fallback for THIS equation only
    }
  }
  
  if (pilot == 'lasso') {
    cv_pilot <- cv.glmnet(X, y_k, alpha=1,
                          nfolds=n_folds, standardize=standardize)
    b_pilot  <- as.numeric(coef(cv_pilot, s='lambda.min'))[-1]
  }
  # Adaptive weights
  w <- 1 / (abs(b_pilot)^nu + weight_eps)
  w <- w / mean(w)

  # Stage 2: adaptive LASSO
  cv_fit <- cv.glmnet(X, y_k,
                      alpha          = alpha,
                      nfolds         = n_folds,
                      nlambda        = nlambda,
                      standardize    = standardize,
                      penalty.factor = w)

  b_sel          <- as.numeric(coef(cv_fit, s=cv_criterion))[-1]
  B[, k]         <- b_sel
  lambda_sel[k]  <- cv_fit[[cv_criterion]]
}

# --- Write output -------------------------------------------------------
con <- file(out_file, 'w')

# Section 1: B matrix
writeLines('B', con)
for (row in 1:nrow(B)) {
  writeLines(paste(B[row,], collapse=','), con)
}

# Section 2: meta
writeLines('meta', con)
writeLines(paste0('lambda,',      paste(lambda_sel,        collapse=',')), con)
writeLines(paste0('selected,',    paste(as.integer(B!=0),  collapse=',')), con)

fitted    <- X %*% B
residuals <- Y0 - fitted
mse       <- colMeans(residuals^2)
writeLines(paste0('mse,',         paste(mse,               collapse=',')), con)
writeLines(paste0('sparsity,',    mean(B==0)),                              con)
writeLines(paste0('n_nonzero,',   sum(B!=0)),                               con)
writeLines(paste0('pilot,',       pilot),                                   con)
writeLines(paste0('cv_criterion,',cv_criterion),                            con)
writeLines(paste0('nu,',          nu),                                      con)
writeLines(paste0('oos_msfe,',    'NaN'),                                   con)

close(con)