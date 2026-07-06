% methods/msvar/msvar_config.m

function cfg = msvar_config()
% Default config for msVAR (modified sparse VAR via TSGlasso) via R.
%
% Reference: Dallakyan, Kim & Pourahmadi (2022), CSDA 176, 107557.
%
% Stage 1: TSGlasso identifies significant pairs from the inverse
%          spectral density matrix (ADMM-based penalised log-likelihood).
% Stage 2: FDR (Benjamini-Hochberg) refines non-zero AR coefficients.

    cfg.method      = 'msvar';
    cfg.p           = [];       % lag order; overridden by BIC selection
    cfg.p_seq       = [1 2 3];  % lag orders to search over (p.seq in R)

    % --- Stage 1: TSGlasso / spectral screening -------------------------
    cfg.lambda              = [];       % [] → auto-selected by selectlambda
    cfg.half_window_length  = 12;       % halfWindowLength (spectral smoother)
    cfg.rho                 = 10;       % ADMM penalty parameter
    cfg.alpha               = 1.5;      % ADMM relaxation
    cfg.thresh              = 1e-6;     % inverse spectral density threshold
    cfg.rho_flex            = true;     % adaptive rho in ADMM
    cfg.ADMM_iter           = 100;      % max ADMM iterations
    cfg.stage1_info         = 'bic';    % 'bic' or 'aic' for lag selection
    cfg.standardize         = false;    % standardize data before fitting
    cfg.stage1_show_status  = false;

    % --- lambda selection (only used when cfg.lambda is empty) ----------
    cfg.lambda_gam          = 0.5;      % eBIC gamma for selectlambda
    cfg.lambda_trim_max     = 0.8;
    cfg.lambda_trim_min     = 0.01;
    cfg.lambda_criteria     = 'BIC';    % 'AIC','BIC','eBIC','CV'

    % --- Stage 2: FDR refinement ----------------------------------------
    cfg.fdr_q               = 0.1;      % FDR significance threshold (q)
    cfg.stage2_show_status  = false;

    % --- Numerical optimisation -----------------------------------------
    cfg.ite_max             = 300;
    cfg.reltol              = 1e-3;

    % --- R interface ----------------------------------------------------
    cfg.r_script_path = '';   % auto-detected if empty
    cfg.r_exe = 'Rscript';
    cfg.tmp_dir       = tempdir();
    cfg.keep_tmp      = false;

    % --- output ---------------------------------------------------------
    cfg.save    = true;
    cfg.verbose = false;
end