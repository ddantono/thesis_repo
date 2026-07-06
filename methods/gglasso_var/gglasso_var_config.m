% methods/gglasso_var/gglasso_var_config.m

function cfg = gglasso_var_config()
% Default config for Group LASSO VAR via R/gglasso.
%
% Groups are defined lag-wise: all K coefficients at a given lag for a
% given equation form one group of size K.  This imposes a structured
% sparsity prior — entire lag-blocks are either included or zeroed.
%
% Reference: Yang & Zou (2015), Statistics and Computing 25(6), 1129-1141.
%            R package: gglasso (CRAN).

    cfg.method      = 'gglasso_var';
    cfg.p           = 1;            % lag order

    % --- group lasso (gglasso) ------------------------------------------
    cfg.nlambda      = 100;         % number of lambda values on path
    cfg.lambda       = [];          % [] → auto-selected by CV
    cfg.loss         = 'ls';        % least squares (regression)
    cfg.cv_criterion = 'lambda.min';% 'lambda.min' or 'lambda.1se'
    cfg.n_folds      = 10;          % k-fold CV folds
    cfg.eps          = 1e-8;        % convergence tolerance
    cfg.maxit        = 3e6;         % max inner-loop iterations
    cfg.intercept    = false;       % intercept per equation (gglasso param)

    % --- grouping scheme ------------------------------------------------
    % 'lag'      : one group per (lag, equation) — K coefficients each
    % 'variable' : one group per (predictor variable) across all lags — p per group
    cfg.group_by     = 'lag';

    % --- R interface ----------------------------------------------------
    cfg.r_script_path = '';         % auto-detected if empty
    cfg.r_exe = 'Rscript';
    cfg.tmp_dir       = tempdir();
    cfg.keep_tmp      = false;

    % --- output ---------------------------------------------------------
    cfg.save    = true;
    cfg.verbose = false;
end