# methods/pds_lm_var/pds_lm_var.R
#
# Post-Double-Selection LM Granger Causality VAR
# Hecq, Margaritella & Smeekes (2021), J. Financial Econometrics
#
# Called from run_pds_lm_var.m via call_r_script().
#
# ARGS (command line):
#   argv[1] : input CSV  — Y matrix (T x K), header y1..yK
#   argv[2] : config CSV — param,value pairs
#   argv[3] : output CSV — sections B / pvalues / gc_matrix / fstats / meta
#
# Requires: glmnet

suppressPackageStartupMessages({
  library(glmnet)
})

# =========================================================================
# 0. Parse arguments
# =========================================================================
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 3) stop("Usage: pds_lm_var.R <input> <config> <output>")

in_file  <- argv[1]
cfg_file <- argv[2]
out_file <- argv[3]

# =========================================================================
# 1. Read inputs
# =========================================================================
Y_raw <- read.csv(in_file, header = TRUE)
Y     <- as.matrix(Y_raw)
T_obs <- nrow(Y)
K     <- ncol(Y)

cfg_raw <- read.csv(cfg_file, header = TRUE, stringsAsFactors = FALSE)
cfg     <- setNames(as.list(cfg_raw$value), cfg_raw$param)

p                        <- as.integer(cfg$p)
crit                     <- tolower(as.character(cfg$crit))
alpha_pen                <- as.numeric(cfg$alpha)
sign_level               <- as.numeric(cfg$sign)
use_fsc                  <- as.integer(cfg$finite_sample_correction) == 1  # finite-sample correction

# =========================================================================
# 2. Build lagged dataset
# =========================================================================
# Returns list: $Y0 [T_eff x K], $X [T_eff x K*p] lag-major
build_lags <- function(Y, p) {
  T_obs <- nrow(Y)
  K     <- ncol(Y)
  T_eff <- T_obs - p
  Y0    <- Y[(p+1):T_obs, , drop = FALSE]
  X     <- matrix(0, T_eff, K * p)
  for (lag in 1:p) {
    cols <- ((lag-1)*K + 1):(lag*K)
    X[, cols] <- Y[(p+1-lag):(T_obs-lag), , drop = FALSE]
  }
  colnames(X) <- paste0("y", rep(1:K, p), "_lag", rep(1:p, each = K))
  list(Y0 = Y0, X = X)
}

lags_out <- build_lags(Y, p)
Y0 <- lags_out$Y0   # [T_eff x K]
X  <- lags_out$X    # [T_eff x K*p]
T_eff <- nrow(Y0)

# =========================================================================
# 3. ic_glmnetboundT
#    Selects lambda via IC with lower bound: at most floor(T_eff/2) vars.
# =========================================================================
ic_glmnetboundT <- function(x, y, crit = "bic", alpha = 1, standardize = FALSE) {
  n <- length(y)
  
  # Fit glmnet path
  fit   <- glmnet(x = x, y = y, alpha = alpha, standardize = standardize)
  coefs <- as.matrix(coef(fit))   # (p+1) x nlambda,  row 1 = intercept
  
  # Apply lower bound: keep only solutions with <=floor(n/2) nonzero vars
  n_nonzero <- colSums(coefs[2:nrow(coefs), , drop = FALSE] != 0)
  keep      <- n_nonzero <= floor(n / 2)
  if (!any(keep)) keep[length(keep)] <- TRUE  # fallback: keep most sparse
  
  coefs  <- coefs[, keep, drop = FALSE]
  lambda <- fit$lambda[keep]
  df     <- colSums(coefs[2:nrow(coefs), , drop = FALSE] != 0)
  
  # Fitted values and residuals for each lambda
  yhat      <- cbind(1, x) %*% coefs
  residuals <- sweep(yhat, 1, y, "-") * (-1)   # y - yhat
  mse       <- colMeans(residuals^2)
  sse       <- colSums(residuals^2)
  
  nvar <- df + 1   # +1 for intercept
  bic  <- n * log(mse) + nvar * log(n)
  aic  <- n * log(mse) + 2 * nvar
  ebic <- n * log(mse) + nvar * log(n) + 2 * nvar * 0.5 * log(ncol(x))
  aicc <- aic + (2 * nvar * (nvar + 1)) / pmax(n - nvar - 1, 1)
  hqc  <- n * log(mse) + 2 * nvar * log(log(n))
  
  ic_vals <- switch(crit,
                    bic  = bic,
                    aic  = aic,
                    ebic = ebic,
                    aicc = aicc,
                    hqc  = hqc,
                    bic  # default
  )
  
  best <- which.min(ic_vals)
  list(
    coefficients = coefs[, best],
    lambda       = lambda[best],
    nvar         = nvar[best],
    selected     = which(coefs[2:nrow(coefs), best] != 0)  # 1-based col index of x
  )
}

# =========================================================================
# 4. PDS-LM test for one pair (cause_idx -> target_idx)
#
#   Returns list: pvalue, fstat, b_coef (post-selection OLS coef vec, length K*p)
# =========================================================================
pds_lm_pair <- function(Y0, X, target_idx, cause_idx, p, K,
                        crit, alpha_pen, sign_level, use_fsc) {
  
  T_eff   <- nrow(Y0)
  y_resp  <- Y0[, target_idx]   # response variable
  
  # Columns of X corresponding to cause variable lags
  # Lag-major layout: col (lag-1)*K + k = y_k at lag
  cause_cols <- ((0:(p-1)) * K) + cause_idx   # length p
  
  # Control columns = all lags EXCEPT cause lags
  control_cols <- setdiff(1:ncol(X), cause_cols)
  X_ctrl  <- X[, control_cols, drop = FALSE]   # [T_eff x (K-1)*p]
  X_cause <- X[, cause_cols,   drop = FALSE]   # [T_eff x p]
  
  # -----------------------------------------------------------------------
  # Step 1: LASSO of y_resp on X_ctrl (IC-tuned, lower bound)
  # -----------------------------------------------------------------------
  lasso_y <- ic_glmnetboundT(X_ctrl, y_resp, crit = crit, alpha = alpha_pen)
  
  # -----------------------------------------------------------------------
  # Step 2: LASSO of each cause lag on X_ctrl
  # -----------------------------------------------------------------------
  sel_cause <- integer(0)
  for (j in 1:p) {
    lasso_xj <- ic_glmnetboundT(X_ctrl, X_cause[, j],
                                crit = crit, alpha = alpha_pen)
    sel_cause <- union(sel_cause, lasso_xj$selected)
  }
  
  # -----------------------------------------------------------------------
  # Step 3: Union of selected variables
  #   Force own lags of target into selection (double-selection safeguard)
  # -----------------------------------------------------------------------
  target_lag_names <- paste0("y", target_idx, "_lag", 1:p)
  own_lag_cols_in_ctrl <- which(colnames(X_ctrl) %in% target_lag_names)
  
  sel_union <- sort(unique(c(lasso_y$selected, sel_cause, own_lag_cols_in_ctrl)))
  
  if (length(sel_union) == 0) {
    # Degenerate: no controls selected — use own lags only
    sel_union <- own_lag_cols_in_ctrl
  }
  
  X_sel <- X_ctrl[, sel_union, drop = FALSE]   # [T_eff x |sel|]
  
  # -----------------------------------------------------------------------
  # Step 4: OLS of y_resp on X_sel → residuals for auxiliary regression
  # -----------------------------------------------------------------------
  tryCatch({
    fit_main <- lm(y_resp ~ X_sel)
    e_hat    <- fit_main$residuals   # [T_eff]
    
    # Auxiliary regression: e_hat ~ X_cause + X_sel  (Step 4 of Algorithm 1)
    X_aux <- cbind(X_cause, X_sel)   # [T_eff x (p + |sel|)]
    fit_aux <- lm(e_hat ~ X_aux)
    R2  <- summary(fit_aux)$r.squared
    
    s_hat <- length(sel_union)   # number of selected controls
    NGC   <- p                   # degrees of freedom = number of cause lags
    
    if (use_fsc) {
      # Finite-sample F-correction (Step 4b)
      df2  <- T_eff - s_hat - NGC   # conservative
      if (df2 <= 0) df2 <- 1
      fstat <- ((T_eff - s_hat - NGC) / NGC) * (R2 / max(1 - R2, 1e-15))
      pval  <- pf(fstat, df1 = NGC, df2 = df2, lower.tail = FALSE)
    } else {
      # Asymptotic chi2 LM (Step 4a)
      fstat <- T_eff * R2
      pval  <- pchisq(fstat, df = NGC, lower.tail = FALSE)
    }
    
    # ------------------------------------------------------------------
    # Post-selection OLS coefficients for the full equation
    # (response ~ all selected controls + cause lags)
    # Returns a vector of length K*p with zeros for non-selected lags
    # ------------------------------------------------------------------
    X_full_sel <- cbind(X_cause, X_sel)
    fit_full   <- lm(y_resp ~ X_full_sel)
    b_full_sel <- coef(fit_full)[-1]  # drop intercept
    
    # Map back to K*p index space
    b_vec <- rep(0, K * p)
    # cause lags → positions cause_cols
    b_vec[cause_cols] <- b_full_sel[1:p]
    # control lags → positions control_cols[sel_union]
    if (length(sel_union) > 0) {
      ctrl_positions <- control_cols[sel_union]
      b_vec[ctrl_positions] <- b_full_sel[(p+1):(p+length(sel_union))]
    }
    
    list(
      pvalue = pval,
      fstat = fstat,
      b_coef = b_vec,
      selected_y = lasso_y$selected,
      selected_cause = sel_cause,
      selected_union = sel_union,
      success = TRUE
    )
    
  }, error = function(e) {
    list(
      pvalue = NA_real_,
      fstat = NA_real_,
      b_coef = rep(0, K * p),
      selected_y = integer(0),
      selected_cause = integer(0),
      selected_union = integer(0),
      success = FALSE
    )
  })
}

# =========================================================================
# 5. Run all K*(K-1) pairs
# =========================================================================
# Result containers
pval_mat  <- matrix(NA_real_, K, K)   # pval_mat[i,j] = p-value j->i
fstat_mat <- matrix(0,        K, K)
gc_mat    <- matrix(0L,       K, K)
B_mat     <- matrix(0,        K*p, K) # B_mat[,i] = coefs for equation i
residuals_mat <- matrix(0,        T_eff, K)

diag(pval_mat) <- -1   # diagonal sentinel -> decoded as NaN in MATLAB
# NA pvalues from failed pairs also written as -2, decoded as NaN in MATLAB

for (target in 1:K) {
  
  y_resp      <- Y0[, target]
  sig_causes  <- integer(0)   # cause indices with significant GC
  sel_controls <- integer(0)  # union of selected control col indices in X
  
  # --- Pass 1: run all K-1 pair tests, collect results and selections ---
  for (cause in 1:K) {
    if (cause == target) next
    
    res <- pds_lm_pair(Y0, X, target, cause, p, K,
                       crit, alpha_pen, sign_level, use_fsc)
    
    pval_mat[target, cause]  <- if (is.finite(res$pvalue))  res$pvalue  else NA_real_
    fstat_mat[target, cause] <- if (is.finite(res$fstat))   res$fstat   else 0
    is_sig <- !is.na(res$pvalue) && res$pvalue < sign_level
    gc_mat[target, cause]    <- if (is_sig) 1L else 0L
    
    if (is_sig) {
      sig_causes <- c(sig_causes, cause)
    }
    
    # Union of control selections from all pairs (for robust refit)
    # res$selected_union indexes into X_ctrl — map back to X columns
    cause_cols_pair  <- ((0:(p-1)) * K) + cause
    control_cols_pair <- setdiff(1:ncol(X), cause_cols_pair)
    if (length(res$selected_union) > 0) {
      sel_controls <- union(sel_controls, control_cols_pair[res$selected_union])
    }
  }
  
  # --- Pass 2: unified OLS refit for this equation ---------------------
  # Regressors: lags of significant causes + selected controls
  # This avoids double-counting controls across multiple cause pairs.
  
  b_eq <- rep(0, K * p)
  
  if (length(sig_causes) == 0) {
    # No significant causes: fit on selected controls only (for residuals)
    # Coefficients remain zero — equation is unpredictable by other series
    if (length(sel_controls) > 0) {
      X_refit <- X[, sel_controls, drop = FALSE]
      fit_refit <- tryCatch(lm(y_resp ~ X_refit),
                            error = function(e) NULL)
      if (!is.null(fit_refit)) {
        b_refit <- coef(fit_refit)[-1]
        b_eq[sel_controls] <- b_refit
      }
    }
  } else {
    # Collect all cause lag columns for significant causes
    sig_cause_cols <- as.vector(sapply(sig_causes,
                                       function(c) ((0:(p-1)) * K) + c))
    # Union with selected controls, remove any overlap
    refit_cols <- sort(unique(c(sig_cause_cols, sel_controls)))
    X_refit    <- X[, refit_cols, drop = FALSE]
    
    fit_refit <- tryCatch(lm(y_resp ~ X_refit), error = function(e) NULL)
    if (!is.null(fit_refit)) {
      b_refit <- coef(fit_refit)[-1]  # drop intercept
      b_eq[refit_cols] <- b_refit
    }
  }
  
  B_mat[, target] <- b_eq
  
  # In-sample residuals from unified refit
  fitted_target <- X %*% matrix(b_eq, ncol = 1)
  residuals_mat[, target] <- y_resp - fitted_target[, 1]
}

# =========================================================================
# 6. Summary meta
# =========================================================================
n_gc_detected <- sum(gc_mat)
n_pairs       <- K * (K - 1)
sparsity_gc   <- 1 - n_gc_detected / n_pairs

# Canonical diagnostics (computed from final B and residuals)
tol        <- 1e-8
selected   <- as.integer(abs(B_mat) > tol)
n_nonzero  <- sum(selected)
sparsity   <- mean(selected == 0)
mse_per_eq <- colMeans(residuals_mat^2)
oos_msfe   <- NaN

# =========================================================================
# 7. Write output CSV
# =========================================================================
# Section format matches read_r_output_file.m:
#   <SECTION_NAME>   ← no comma
#   rows of comma-separated numbers
#
# Section B       : K*p rows x K cols
# Section pvalues : K rows x K cols  (diagonal = -1 as NaN sentinel)
# Section gc_matrix: K rows x K cols
# Section fstats  : K rows x K cols
# Section meta    : key,value
# Mandatory diagnostics: sparsity, n_nonzero, mse_per_eq, oos_msfe
# Replace R NA with -2 sentinel before writing (MATLAB reads as numeric, decodes to NaN)
pval_write  <- pval_mat
fstat_write <- fstat_mat
pval_write[is.na(pval_write)]   <- -2
fstat_write[is.na(fstat_write)] <- 0

fid <- file(out_file, "w")

# --- B ------------------------------------------------------------------
cat("B\n", file = fid)
for (r in 1:nrow(B_mat)) {
  cat(paste(B_mat[r, ], collapse = ","), "\n", file = fid, sep = "")
}

# --- pvalues ------------------------------------------------------------
cat("pvalues\n", file = fid)
for (r in 1:K) {
  cat(paste(pval_write[r, ], collapse = ","), "\n", file = fid, sep = "")
}

# --- gc_matrix ----------------------------------------------------------
cat("gc_matrix\n", file = fid)
for (r in 1:K) {
  cat(paste(gc_mat[r, ], collapse = ","), "\n", file = fid, sep = "")
}

# --- fstats -------------------------------------------------------------
cat("fstats\n", file = fid)
for (r in 1:K) {
  cat(paste(fstat_write[r, ], collapse = ","), "\n", file = fid, sep = "")
}

# --- meta ---------------------------------------------------------------
cat("meta\n", file = fid)

# Canonical diagnostics schema (mandatory fields first)
cat(sprintf("sparsity,%.10g\n",    sparsity),                                file = fid)
cat(sprintf("n_nonzero,%d\n",      as.integer(n_nonzero)),                   file = fid)
cat(sprintf("mse_per_eq,%s\n",     paste(sprintf("%.15g", mse_per_eq), collapse = ",")), file = fid)
cat(sprintf("oos_msfe,NaN\n"),                                              file = fid)

# Method-specific metadata
cat(sprintf("n_gc_detected,%d\n",  as.integer(n_gc_detected)),               file = fid)
cat(sprintf("n_pairs,%d\n",        as.integer(n_pairs)),                     file = fid)
cat(sprintf("sparsity_gc,%.10g\n", sparsity_gc),                             file = fid)
cat(sprintf("selected,%s\n",       paste(selected, collapse = ",")),         file = fid)
cat(sprintf("K,%d\n",              as.integer(K)),                           file = fid)
cat(sprintf("p,%d\n",              as.integer(p)),                           file = fid)
cat(sprintf("T_eff,%d\n",          as.integer(T_eff)),                       file = fid)
cat(sprintf("sign,%.10g\n",        sign_level),                              file = fid)
cat(sprintf("crit,%s\n",           crit),                                    file = fid)
close(fid)

cat(sprintf("[pds_lm_var] Done. K=%d, p=%d, GC pairs detected: %d/%d\n",
            K, p, n_gc_detected, n_pairs))