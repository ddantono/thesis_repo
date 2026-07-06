function cfg = lasso_var_config()
% Default config for LASSO-VAR and Elastic Net-VAR.
% Set cfg.alpha = 1 for pure LASSO, 0 < alpha < 1 for Elastic Net.

    cfg.method      = 'lasso_enet_var';
    cfg.p           = 1;          % VAR lag order

    % --- regularisation --------------------------------------------------
    cfg.lambda      = [];         % [] → CV selects lambda automatically
    cfg.alpha       = 1.0;        % 1 = LASSO, (0,1) = Elastic Net

    % --- cross-validation ------------------------------------------------
    cfg.cv          = true;       % use CV to select lambda
    cfg.n_folds     = 10;
    cfg.cv_criterion = 'MSE';     % 'MSE' or 'min+1se'

    % --- estimation ------------------------------------------------------
    cfg.standardize = true;       % standardize regressors inside lasso()
    cfg.max_iter    = 1e4;
    cfg.rel_tol     = 1e-4;

    % --- output ----------------------------------------------------------
    cfg.save        = true;
    cfg.verbose     = false;
end