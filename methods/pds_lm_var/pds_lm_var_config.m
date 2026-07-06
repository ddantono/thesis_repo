% methods/pds_lm_var/pds_lm_var_config.m

function cfg = pds_lm_var_config()
% Default config for Post-Double-Selection LM Granger Causality VAR
% via R (Hecq, Margaritella & Smeekes, 2021, J. Financial Econometrics).
%
% Tests Granger causality for all K*(K-1) directed pairs in a K-dim VAR.
% Uses BIC-tuned LASSO with lower bound on penalty (max floor(T/2) vars).
%
% Reference: Hecq, Margaritella & Smeekes (2021),
%   "Granger Causality Testing in High-Dimensional VARs:
%    A Post-Double-Selection Procedure", J. Financial Econometrics.

    cfg.method  = 'pds_lm_var';
    cfg.p       = 1;              % VAR lag order

    % --- penalization ---------------------------------------------------
    cfg.crit        = 'bic';     % IC for lambda: 'bic','aic','ebic','aicc','hqc'
    cfg.alpha       = 1.0;       % 1=LASSO, 0.5=Elastic Net (glmnet alpha)
    cfg.standardize = false;     % standardize regressors inside ic_glmnetboundT

    % --- test -----------------------------------------------------------
    cfg.sign        = 0.05;      % significance level for GC test
    cfg.finite_sample_correction = true;  % use F-correction (Step 4b) vs chi2 (4a)

    % --- R interface ----------------------------------------------------
    cfg.r_script_path = '';      % auto-detected if empty
    cfg.r_exe = 'Rscript';
    cfg.tmp_dir       = tempdir();
    cfg.keep_tmp      = false;

    % --- output ---------------------------------------------------------
    cfg.save    = true;
    cfg.verbose = false;
end