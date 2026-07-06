% methods/adaptive_lasso_var/adaptive_lasso_var_config.m

function cfg = adaptive_lasso_var_config()
% Default config for Adaptive LASSO-VAR via R/glmnet.

    cfg.method      = 'adaptive_lasso_var';
    cfg.p           = 1;

    % --- pilot estimation -----------------------------------------------
    cfg.pilot       = 'ols';      % 'ols' or 'lasso'
    cfg.nu          = 1;          % weight exponent: w = 1/(|beta|^nu + eps)
    cfg.weight_eps  = 1e-6;

    % --- adaptive LASSO (second stage) ----------------------------------
    cfg.cv_criterion = 'lambda.min';  % 'lambda.min' or 'lambda.1se'
    cfg.n_folds      = 10;
    cfg.nlambda      = 100;
    cfg.standardize  = true;
    cfg.alpha        = 1.0;

    % --- R interface ----------------------------------------------------
    cfg.r_script_path = '';       % auto-detected if empty
    cfg.r_exe = 'Rscript';
    cfg.tmp_dir       = tempdir();% directory for I/O files
    cfg.keep_tmp      = false;    % keep temp files after run

    % --- output ---------------------------------------------------------
    cfg.save    = true;
    cfg.verbose = false;
end