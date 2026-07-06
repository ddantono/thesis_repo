% methods/bgr_var/bgr_var_config.m

function cfg = bgr_var_config()
% Default config for Bayesian Ridge Regression VAR (BGR) via R/BigVAR.
%
%   Implements the Bayesian Ridge (Minnesota-style) VAR from BigVAR.
%   struct = 'BGR' in BigVAR — Banbura, Giannone & Reichlin (2010).
%
%   Reference: Banbura, Giannone & Reichlin, J. Applied Econometrics 2010.
%              Nicholson, Matteson & Bien (2017), BigVAR R package.

    cfg.method  = 'bgr_var';
    cfg.p       = 1;              % VAR lag order

    % --- BigVAR penalty grid --------------------------------------------
    cfg.nlambda = 50;             % grid depth  (gran[1] in BigVAR)
    cfg.n_folds = 10;             % rolling CV evaluations (gran[2])

    % --- BigVAR CV split ------------------------------------------------
    % [] → auto-computed in wrapper as floor(T/3), floor(2*T/3)
    cfg.T1      = [];
    cfg.T2      = [];

    % --- Minnesota prior ------------------------------------------------
    % BGR uses an implicit Minnesota-style shrinkage.
    % Set Minnesota=TRUE to add explicit prior weighting.
    cfg.Minnesota = false;

    % --- intercept ------------------------------------------------------
    cfg.intercept = true;

    % --- R interface ----------------------------------------------------
    cfg.r_script_path = '';       % auto-detected if empty
    cfg.r_exe = 'Rscript';
    cfg.tmp_dir       = tempdir();
    cfg.keep_tmp      = false;

    % --- output ---------------------------------------------------------
    cfg.save    = true;
    cfg.verbose = false;
end
