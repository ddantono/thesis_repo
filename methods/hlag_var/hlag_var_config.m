% methods/hlag_var/hlag_var_config.m

function cfg = hlag_var_config()
% Default config for HLAG-VAR via R/BigVAR.
%
%   Supports struct = 'HLAGOO' (Own/Other) or 'HLAGC' (Componentwise).
%   Both share the same wrapper run_hlag_var.m — only cfg.struct differs.
%
%   Reference: Nicholson, Wilms, Bien & Matteson (JMLR 2020).

    cfg.method  = 'hlag_var';      % dispatch key
    cfg.p       = 1;               % maximum lag order

    % --- BigVAR penalty structure ---------------------------------------
    cfg.struct  = 'HLAGOO';        % 'HLAGOO' or 'HLAGC'

    % --- penalty grid ---------------------------------------------------
    cfg.nlambda = 50;              % grid depth  (gran[1] in BigVAR)
    cfg.n_folds = 10;              % rolling CV evaluations (gran[2])

    % --- BigVAR CV split ------------------------------------------------
    % T1, T2: indices for CV start / forecast eval start.
    % [] → auto-computed in wrapper as floor(T/3), floor(2*T/3).
    cfg.T1      = [];
    cfg.T2      = [];

    % --- relaxed VAR (RVAR) ---------------------------------------------
    % RVAR=TRUE refits OLS on the selected support → cleaner sparsity.
    cfg.RVAR    = true;

    % --- Minnesota prior ------------------------------------------------
    cfg.Minnesota = false;

    % --- intercept ------------------------------------------------------
    cfg.intercept = true;

    % hlag_var_config.m — πρόσθεσε:
    cfg.rel_threshold = 0;   % 0 = καμία post-processing (raw BigVAR output)
                            % 0.05 = 5% of max(|B|) για sparsity visualization

    % --- R interface ----------------------------------------------------
    cfg.r_script_path = '';        % auto-detected if empty
    cfg.r_exe = 'Rscript';
    cfg.tmp_dir       = tempdir();
    cfg.keep_tmp      = false;

    % --- output ---------------------------------------------------------
    cfg.save    = true;
    cfg.verbose = false;
end