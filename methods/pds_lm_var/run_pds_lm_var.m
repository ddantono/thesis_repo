% methods/pds_lm_var/run_pds_lm_var.m

function out = run_pds_lm_var(Y, cfg)
% RUN_PDS_LM_VAR  Post-Double-Selection LM Granger Causality VAR via R.
%
%   Runs the PDS-LM test of Hecq, Margaritella & Smeekes (2021) for all
%   K*(K-1) directed pairs in a K-dimensional VAR(p).
%
%   Each pair (x->y) is tested by:
%     1. LASSO of y on controls (all lags except x lags) — IC-tuned
%     2. LASSO of each x_lag on controls — IC-tuned
%     3. Union of selected variables + forced own lags
%     4. OLS on union → residuals
%     5. Auxiliary regression → F-statistic (finite-sample) or LM (chi2)
%
%   Requires: R on PATH, glmnet R package.
%   Reference: Hecq, Margaritella & Smeekes (2021), J. Financial Econometrics.
%
%   INPUT
%     Y   : [T x K] double — stationary multivariate time series
%     cfg : struct from pds_lm_var_config()
%
%   OUTPUT
%     out : canonical struct — key fields:
%       .coefficients         [K*p x K]  post-selection OLS coef matrix
%       .diagnostics.gc_matrix [K x K]   binary Granger causality matrix
%                                         gc_matrix(i,j)=1 → yj GC-causes yi
%       .diagnostics.pvalue_matrix [K x K]  p-values (NaN on diagonal)
%       .diagnostics.fstat_matrix  [K x K]  F/LM statistics
%       .hyperparameters.crit       char
%       .hyperparameters.alpha      scalar
%       .hyperparameters.sign       scalar

    t_start = tic;
    METHOD  = 'pds_lm_var';

    % ------------------------------------------------------------------ %
    %  1. Config                                                           %
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(cfg)
        cfg = pds_lm_var_config();
    end
    cfg = i_fill_defaults(cfg);

    % ------------------------------------------------------------------ %
    %  2. Input guards                                                     %
    % ------------------------------------------------------------------ %
    if ~isnumeric(Y) || ~ismatrix(Y) || ~isreal(Y)
        out = i_early_failure(METHOD, cfg, t_start, ...
            'Y must be a real 2D numeric matrix.');
        return
    end
    if any(~isfinite(Y(:)))
        out = i_early_failure(METHOD, cfg, t_start, ...
            'Y contains NaN or Inf values.');
        return
    end

    [T, K] = size(Y);
    if isfield(cfg, 'pmax')
        p = cfg.pmax;
    elseif isfield(cfg, 'p')
        p = cfg.p;
    else
        error('run_lasso_enet_var: no lag order found in cfg');
    end

    if K < 2
        out = i_early_failure(METHOD, cfg, t_start, ...
            'PDS-LM requires K>=2 variables for Granger causality testing.');
        return
    end
    if T <= p + 1
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Insufficient observations T=%d for p=%d.', T, p));
        return
    end
    % Need enough obs for post-selection OLS: T_eff > K*p (conservative)
    T_eff = T - p;
    if T_eff <= K * p
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('T_eff=%d <= K*p=%d: too few observations for PDS-LM.', ...
                    T_eff, K*p));
        return
    end

    % ------------------------------------------------------------------ %
    %  3. Check R                                                          %
    % ------------------------------------------------------------------ %
    [r_ok, r_msg] = i_check_r(cfg.r_exe);
    if ~r_ok
        out = i_early_failure(METHOD, cfg, t_start, r_msg);
        return
    end

    % ------------------------------------------------------------------ %
    %  4. Locate R script                                                  %
    % ------------------------------------------------------------------ %
    r_script = i_resolve_r_script(cfg);
    if isempty(r_script)
        out = i_early_failure(METHOD, cfg, t_start, ...
            'pds_lm_var.R not found. Check cfg.r_script_path.');
        return
    end

    % ------------------------------------------------------------------ %
    %  5. Temp file paths                                                  %
    % ------------------------------------------------------------------ %
    run_id   = sprintf('%d', round(now * 1e6));
    in_file  = fullfile(cfg.tmp_dir, sprintf('pdslm_input_%s.csv',  run_id));
    cfg_file = fullfile(cfg.tmp_dir, sprintf('pdslm_config_%s.csv', run_id));
    out_file = fullfile(cfg.tmp_dir, sprintf('pdslm_output_%s.csv', run_id));

    % ------------------------------------------------------------------ %
    %  6. Write input files                                                %
    % ------------------------------------------------------------------ %
    try
        write_r_input_file(Y, in_file);

        cfg_r                          = struct();
        cfg_r.p                        = p;
        cfg_r.K                        = K;
        cfg_r.T                        = T;
        cfg_r.crit                     = cfg.crit;
        cfg_r.alpha                    = cfg.alpha;
        cfg_r.standardize              = cfg.standardize;
        cfg_r.sign                     = cfg.sign;
        cfg_r.finite_sample_correction = cfg.finite_sample_correction;

        write_r_config_file(cfg_r, cfg_file);

    catch ME
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Failed to write input files: %s', ME.message));
        return
    end

    % ------------------------------------------------------------------ %
    %  7. Call R                                                           %
    % ------------------------------------------------------------------ %
    [r_ok, r_log] = call_r_script(cfg.r_exe, r_script, ...
                                   in_file, cfg_file, out_file);
    if ~r_ok
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('R call failed: %s', r_log));
        i_cleanup(cfg, in_file, cfg_file, out_file);
        return
    end

    % ------------------------------------------------------------------ %
    %  8. Read R output                                                    %
    % ------------------------------------------------------------------ %
    try
        r_result = read_r_output_file(out_file);
    catch ME
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Failed to read R output: %s', ME.message));
        i_cleanup(cfg, in_file, cfg_file, out_file);
        return
    end

    % ------------------------------------------------------------------ %
    %  9. Cleanup                                                          %
    % ------------------------------------------------------------------ %
    i_cleanup(cfg, in_file, cfg_file, out_file);

    % ------------------------------------------------------------------ %
    %  10. Validate R output sections                                      %
    % ------------------------------------------------------------------ %
    [r_ok, r_msg] = i_validate_r_result(r_result, K, p);
    if ~r_ok
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('R output validation failed: %s', r_msg));
        return
    end

    % ------------------------------------------------------------------ %
    %  11. Normalize → canonical output                                   %
    % ------------------------------------------------------------------ %
    out = i_normalize_pds_lm(r_result, METHOD, cfg, t_start, T, K, p, Y);

    % ------------------------------------------------------------------ %
    %  12. Post-validation                                                 %
    % ------------------------------------------------------------------ %
    out = i_pds_lm_postvalidate(out, K, p);
end


% ======================================================================= %
%  Private: normalization                                                  %
% ======================================================================= %

function out = i_normalize_pds_lm(r_result, method_name, cfg, t_start, T, K, p, Y)
% Converts PDS-LM R output to canonical struct.
%
% R output sections expected:
%   B          : [K*p x K]   post-selection OLS coefficient matrix
%                             assembled column-by-column (eq i = col i)
%   pvalues    : [K x K]     p-value matrix (row=target, col=cause)
%                             diagonal = NaN encoded as -1 in CSV
%   gc_matrix  : [K x K]     binary GC matrix (same orientation)
%   fstats     : [K x K]     F/LM statistics (diagonal = 0)
%   meta       : key,value   scalar summary fields

    T_eff = T - p;

    % --- Extract B ------------------------------------------------------
    B = r_result.B;   % [K*p x K]

    % --- Fitted values and residuals ------------------------------------
    [Y0, X] = i_build_lagged(Y, T, K, p);
    fitted    = X * B;
    residuals = Y0 - fitted;

    warn_msgs = {};
    if any(~isfinite(residuals(:)))
        warn_msgs{end+1} = 'Residuals contain NaN or Inf.';
    end

    % --- GC results matrices --------------------------------------------
    pval_mat = r_result.pvalues;   % [K x K]
    gc_mat   = r_result.gc_matrix; % [K x K]
    fstat_mat= r_result.fstats;    % [K x K]

    % Decode diagonal sentinel (-1 → NaN)
    pval_mat(pval_mat < 0) = NaN;

    % --- Meta fields ----------------------------------------------------
    meta_r = i_safe_struct(r_result, 'meta');

    n_gc_detected = i_meta_num(meta_r, 'n_gc_detected', NaN);
    n_pairs       = i_meta_num(meta_r, 'n_pairs',       K*(K-1));
    sparsity_gc   = i_meta_num(meta_r, 'sparsity_gc',   NaN);

    % --- Sparsity of B --------------------------------------------------
    sparsity_B = sum(B(:) == 0) / numel(B);
    n_nonzero  = sum(B(:) ~= 0);
    mse_per_eq = mean(residuals.^2, 1);   % [1 x K]

    % --- Pack raw struct -------------------------------------------------
    raw                  = struct();
    raw.B                = B;
    raw.Y0               = Y0;
    raw.residuals        = residuals;
    raw.fitted_values    = fitted;
    raw.T                = T;
    raw.K                = K;
    raw.T_eff            = T_eff;
    raw.p                = p;
    raw.lambda           = NaN;   % no single lambda — IC-selected per regression
    raw.toolbox_used     = 'R:glmnet';
    raw.selected_terms   = i_selected_term_labels(B, K, p);
    raw.warnings         = warn_msgs;

    raw.diagnostics = struct( ...
        'sparsity',        sparsity_B, ...
        'n_nonzero',       n_nonzero, ...
        'mse_per_eq',      mse_per_eq, ...
        'oos_msfe',        NaN, ...
        'gc_matrix',       gc_mat, ...
        'pvalue_matrix',   pval_mat, ...
        'fstat_matrix',    fstat_mat, ...
        'n_gc_detected',   n_gc_detected, ...
        'n_pairs_tested',  n_pairs, ...
        'sparsity_gc',     sparsity_gc);

    raw.hyperparameters = struct( ...
        'crit',                     cfg.crit, ...
        'alpha',                    cfg.alpha, ...
        'sign',                     cfg.sign, ...
        'standardize',              cfg.standardize, ...
        'finite_sample_correction', cfg.finite_sample_correction);

    out = normalize_method_output(raw, method_name, cfg, t_start);
end


% ======================================================================= %
%  Private: validation                                                     %
% ======================================================================= %

function [ok, msg] = i_validate_r_result(r_result, K, p)
    ok  = true;
    msg = '';

    required = {'B', 'pvalues', 'gc_matrix', 'fstats'};
    for i = 1:numel(required)
        if ~isfield(r_result, required{i})
            ok  = false;
            msg = sprintf('Missing R output section: "%s".', required{i});
            return
        end
    end

    % B dimensions
    if ~isequal(size(r_result.B), [K*p, K])
        ok  = false;
        msg = sprintf('B size [%d x %d] expected [%d x %d].', ...
                      size(r_result.B,1), size(r_result.B,2), K*p, K);
        return
    end

    % pvalues, gc_matrix, fstats must be [K x K]
    for fn = {'pvalues', 'gc_matrix', 'fstats'}
        f = fn{1};
        if ~isequal(size(r_result.(f)), [K, K])
            ok  = false;
            msg = sprintf('%s size [%d x %d] expected [%d x %d].', ...
                          f, size(r_result.(f),1), size(r_result.(f),2), K, K);
            return
        end
    end
end

function out = i_pds_lm_postvalidate(out, K, p)
    if ~out.success; return; end

    % Warn if GC matrix missing
    if ~isfield(out.diagnostics, 'gc_matrix') || ...
            isempty(out.diagnostics.gc_matrix)
        out.warnings{end+1} = 'gc_matrix is empty — no GC results.';
        out.status_code = 2;
    end

    % Warn if all p-values are NaN
    if isfield(out.diagnostics, 'pvalue_matrix')
        pv = out.diagnostics.pvalue_matrix;
        off_diag = pv(~isnan(pv));
        if isempty(off_diag) || all(isnan(off_diag))
            out.warnings{end+1} = 'All p-values are NaN.';
            out.status_code = 2;
        end
    end

    % Warn if all B coefficients zero (post-selection OLS degenerate)
    if ~isempty(out.coefficients) && all(out.coefficients(:) == 0)
        out.warnings{end+1} = ...
            'All post-selection OLS coefficients zero.';
        out.status_code = 2;
    end

    % Document that PDS-LM B is post-selection OLS — not penalized estimation
    out.warnings{end+1} = ...
        'PDS-LM B is post-selection OLS refit — use gc_matrix for GC evaluation, not B/mse/sparsity.';
end


% ======================================================================= %
%  Private: helpers                                                        %
% ======================================================================= %

function [Y0, X] = i_build_lagged(Y, T, K, p)
    Y0 = Y(p+1:end, :);
    X  = zeros(T-p, K*p);
    for lag = 1:p
        cols       = (lag-1)*K + (1:K);
        X(:, cols) = Y(p+1-lag:end-lag, :);
    end
end

function labels = i_selected_term_labels(B, K, p)
    nz     = any(B ~= 0, 2);
    idx    = find(nz);
    labels = cell(numel(idx), 1);
    for i = 1:numel(idx)
        lag_num   = ceil(idx(i) / K);
        var_num   = mod(idx(i)-1, K) + 1;
        labels{i} = sprintf('y%d_lag%d', var_num, lag_num);
    end
end

function v = i_meta_num(meta, fname, default)
    if isstruct(meta) && isfield(meta, fname) && ~isempty(meta.(fname))
        v = meta.(fname);
        if isnumeric(v); v = v(1); end
    else
        v = default;
    end
end

function s = i_safe_struct(r_result, fname)
    if isfield(r_result, fname) && isstruct(r_result.(fname))
        s = r_result.(fname);
    else
        s = struct();
    end
end

function cfg = i_fill_defaults(cfg)
    defaults      = pds_lm_var_config();
    caller_method = cfg.method;
    fields        = fieldnames(defaults);
    for i = 1:numel(fields)
        f = fields{i};
        if ~isfield(cfg, f)
            cfg.(f) = defaults.(f);
        end
    end
    cfg.method = caller_method;
    if isempty(cfg.p); cfg.p = 1; end
end

function [ok, msg] = i_check_r(r_exe)
    [status, ~] = system(sprintf('"%s" --version', r_exe));
    if status == 0
        ok  = true;
        msg = '';
    else
        ok  = false;
        msg = sprintf('R not found at "%s".', r_exe);
    end
end

function r_script = i_resolve_r_script(cfg)
    if ~isempty(cfg.r_script_path) && exist(cfg.r_script_path, 'file')
        r_script = cfg.r_script_path;
        return
    end
    this_dir  = fileparts(mfilename('fullpath'));
    candidate = fullfile(this_dir, 'pds_lm_var.R');
    if exist(candidate, 'file')
        r_script = candidate;
    else
        r_script = '';
    end
end

function i_cleanup(cfg, varargin)
    if cfg.keep_tmp; return; end
    for i = 1:numel(varargin)
        f = varargin{i};
        if exist(f, 'file'), delete(f); end
    end
end

function out = i_early_failure(method, cfg, t_start, msg)
    raw         = struct();
    raw.success = false;
    raw.message = msg;
    out         = normalize_method_output(raw, method, cfg, t_start);
end