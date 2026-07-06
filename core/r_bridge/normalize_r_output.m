% core/r_bridge/normalize_r_output.m

function out = normalize_r_output(r_result, method_name, cfg, ...
                                   t_start, T, K, p, Y)
% NORMALIZE_R_OUTPUT  Converts any R result struct to canonical output.
%
%   r_result : struct from read_r_output_file — must have .B and .meta
%   Y        : [T x K] original data (optional — pass [] if unavailable)
%
%   Fitted values and residuals are recomputed from B and Y if Y provided.
%   All missing meta fields default to sentinel values.

    T_eff = T - p;

    % --- Extract B ------------------------------------------------------
    if ~isfield(r_result, 'B') || isempty(r_result.B)
        raw         = struct('success', false, ...
                             'message', 'R output missing B matrix.');
        out         = normalize_method_output(raw, method_name, cfg, t_start);
        return
    end

    B    = r_result.B;
    meta = i_safe_struct(r_result, 'meta');

    % --- Sanity check B dimensions --------------------------------------
    warn_msgs = {};
    if ~isequal(size(B), [K*p, K])
        warn_msgs{end+1} = sprintf( ...
            'B size [%d x %d] does not match expected [%d x %d].', ...
            size(B,1), size(B,2), K*p, K);
    end

    % --- Fitted values and residuals ------------------------------------
    if ~isempty(Y) && isnumeric(Y)
        [Y0, X] = i_build_lagged(Y, T, K, p);
        fitted    = X * B;
        residuals = Y0 - fitted;

        % NaN/Inf check
        if any(~isfinite(residuals(:)))
            warn_msgs{end+1} = 'Residuals contain NaN or Inf.';
        end
    else
        Y0        = [];
        fitted    = [];
        residuals = [];
    end

     % --- Recompute diagnostics from MATLAB B/residuals ------------------
    if ~isempty(residuals)
        matlab_mse_per_eq = mean(residuals.^2, 1);
    else
        matlab_mse_per_eq = i_meta_field(meta, 'mse_per_eq', ...
                                i_meta_field(meta, 'mse', zeros(1,K)));
    end
    matlab_sparsity  = sum(B(:)==0) / numel(B);
    matlab_n_nonzero = sum(B(:)~=0);

    % --- Pack raw struct -------------------------------------------------
    raw                  = struct();
    raw.B                = B;
    raw.Y0               = [];
    raw.residuals        = residuals;
    raw.fitted_values    = fitted;
    raw.T                = T;
    raw.K                = K;
    raw.T_eff            = T_eff;
    raw.p                = p;
    raw.lambda           = i_meta_field(meta, 'lambda',    NaN);
    raw.toolbox_used     = i_cfg_field(cfg, 'toolbox_used', 'R');
    raw.selected_terms   = i_selected_term_labels(B, K, p);
    raw.warnings         = warn_msgs;

    % diagnostics — pull from meta, default to sentinels
    % normalize_r_output.m — στο raw.diagnostics πρόσθεσε:
   raw.diagnostics = struct( ...
    'sparsity',      matlab_sparsity, ...
    'n_nonzero',     matlab_n_nonzero, ...
    'mse_per_eq',    matlab_mse_per_eq, ...
    'lambda_per_eq', i_meta_field(meta, 'lambda',     NaN), ...
    'oos_msfe',      i_meta_field(meta, 'oos_msfe',   NaN), ...
    'pilot',         i_meta_field(meta, 'pilot',      ''));

    % hyperparameters — pull from meta + cfg
    raw.hyperparameters = struct( ...
        'lambda',       i_meta_field(meta, 'lambda',       NaN), ...
        'cv_criterion', i_meta_field(meta, 'cv_criterion', ''), ...
        'n_folds',      i_cfg_field(cfg,  'n_folds',       NaN), ...
        'standardize',  i_cfg_field(cfg,  'standardize',   true), ...
        'alpha',        i_cfg_field(cfg,  'alpha',         1.0));

    % Allow method-specific extra hyperparameters from meta
    extra_keys = setdiff(fieldnames(meta), ...
    {'lambda','mse','mse_per_eq','sparsity','n_nonzero','cv_criterion','selected', ...
     'oos_msfe','pilot'});
    for i = 1:numel(extra_keys)
        k = extra_keys{i};
        raw.hyperparameters.(k) = meta.(k);
    end

    out = normalize_method_output(raw, method_name, cfg, t_start);
end


% --- private helpers ----------------------------------------------------

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

function v = i_meta_field(meta, fname, default)
    if isstruct(meta) && isfield(meta, fname) && ~isempty(meta.(fname))
        v = meta.(fname);
    else
        v = default;
    end
end

function v = i_cfg_field(cfg, fname, default)
    if isfield(cfg, fname) && ~isempty(cfg.(fname))
        v = cfg.(fname);
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