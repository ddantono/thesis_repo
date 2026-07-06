% methods/ncpen_var/ncpen_var_config.m

function cfg = ncpen_var_config()
% Default config for Nonconcave Penalized VAR (SCAD/MCP) via R/bigVAR.
%
% Reference: Davis, Zang & Zheng (2016), Sparse Vector Autoregressive Modeling,
%            Journal of Computational and Graphical Statistics.
%            R package: bigVAR (Nicholson, Matteson & Bien, 2017)

    cfg.method      = 'ncpen_var';
    cfg.p           = 1;          % VAR lag order

    % --- penalty --------------------------------------------------------
    cfg.penalty     = 'SCAD';     % 'SCAD' or 'MCP'
    cfg.lambda      = [];         % [] = CV selects automatically
    cfg.n_folds     = 10;
    cfg.nlambda     = 50;

    % --- bigVAR specific ------------------------------------------------
    cfg.Minnesota   = false;      % Minnesota prior weighting
    cfg.verbose     = false;

    % --- R interface ----------------------------------------------------
    cfg.r_script_path = '';
    cfg.r_exe         = 'Rscript';
    cfg.tmp_dir       = tempdir();
    cfg.keep_tmp      = false;

    % --- output ---------------------------------------------------------
    cfg.save        = true;
end