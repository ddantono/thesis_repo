function out = run_mbts_var(Y, cfg)
% RUN_MBTS_VAR  mBTS sparse VAR wrapper.
%
%   Estimates a sparse VAR model equation-by-equation using the modified
%   Backward-in-Time Selection (mBTS) algorithm (Siggiridou & Kugiumtzis,
%   IEEE TSP 2016). Calls mBTS.m directly — no reimplementation.
%
%   REQUIRES: mBTS.m and DRfitmse.m on the MATLAB path.
%
%   INPUT
%     Y   : [T x K] double — stationary multivariate time series
%     cfg : struct from mbts_var_config() or user-supplied
%
%   OUTPUT
%     out : canonical struct (see create_canonical_output)

    t_start = tic;
    METHOD  = 'mbts_var';

    % ------------------------------------------------------------------ %
    %  1. Config fallback                                                  %
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(cfg)
        cfg = mbts_var_config();
    end
    cfg = i_fill_defaults(cfg);

    % ------------------------------------------------------------------ %
    % ------------------------------------------------------------------ %
    %  2. Dependency check                                                %
    % ------------------------------------------------------------------ %
    if ~exist('DRfitmse', 'file')
        out = i_early_failure(METHOD, cfg, t_start, ...
            'DRfitmse.m not found on MATLAB path. Add methods/mbts_var to path.');
        return
    end
    if ~exist('multilagmatrix', 'file')
        out = i_early_failure(METHOD, cfg, t_start, ...
            'multilagmatrix.m not found on MATLAB path. Add methods/mbts_var to path.');
        return
    end
    % ------------------------------------------------------------------ %
    %  3. Input guards                                                     %
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
    pmax   = cfg.pmax;

    if T <= pmax + 1
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Insufficient observations T=%d for pmax=%d.', T, pmax));
        return
    end
    if pmax < 1 || floor(pmax) ~= pmax
        out = i_early_failure(METHOD, cfg, t_start, ...
            'cfg.pmax must be a positive integer.');
        return
    end

    % ------------------------------------------------------------------ %
    %  4. Run mBTS equation-by-equation                                   %
    %     mBTS(xM, responseindex, pmax) → [indexV, maxorder, MSEval]     %
    %     indexV : [1 x K*pmax] binary — selected lag terms               %
    % ------------------------------------------------------------------ %
    indexM    = zeros(K, K*pmax);   % one row per equation
    maxorders = zeros(1, K);
    mse_per_eq= zeros(1, K);
    warn_msgs = {};

    for k = 1:K
        try
            [indexV, maxorder, MSEval] = mBTS(Y, k, pmax);
            indexM(k, :)   = indexV;
            maxorders(k)   = maxorder;
            mse_per_eq(k)  = MSEval;
        catch ME
            warn_msgs{end+1} = sprintf( ...
                'mBTS failed for equation %d: %s', k, ME.message); %#ok<AGROW>
            % row k of indexM stays zero
        end
    end

    % ------------------------------------------------------------------ %
    %  5. Reconstruct B matrix from indexM                               %
    %     Framework convention: B is [K*pmax x K]                        %
    %     B(lag_var_idx, eq) — non-zero where mBTS selected the term     %
    %     indexM(eq, :) is [1 x K*pmax] in framework lag-major order     %
    % ------------------------------------------------------------------ %
    % mBTS already reshapes indexV to framework-compatible lag-major order
    % (all vars at lag 1, then lag 2, etc.) before returning.
    % So indexM(k,:) maps directly to column k of B.
    B = indexM';   % [K*pmax x K] — binary selection mask

    % ------------------------------------------------------------------ %
    %  6. Compute fitted values and residuals using selected terms        %
    % ------------------------------------------------------------------ %
    p_max_selected = max(maxorders);
    if p_max_selected == 0
        % No terms selected in any equation — flag but continue
        warn_msgs{end+1} = 'mBTS selected no terms in any equation.';
        p_max_selected = 1;
    end

    [Y0, X_full] = i_build_lagged(Y, T, K, pmax);

    % OLS refit on selected terms per equation for coefficients and residuals
    B_coef    = zeros(K*pmax, K);
    fitted    = zeros(T - pmax, K);
    residuals = zeros(T - pmax, K);

    for k = 1:K
        sel = logical(B(:, k));   % selected regressors for equation k
        if any(sel)
            X_sel      = X_full(:, sel);
            b_sel      = X_sel \ Y0(:, k);       % OLS on selected terms
            B_coef(sel, k)   = b_sel;
            fitted(:, k)     = X_sel * b_sel;
        else
            % No regressors selected: fitted = mean(Y0(:,k))
            fitted(:, k) = mean(Y0(:, k));
        end
        residuals(:, k) = Y0(:, k) - fitted(:, k);
    end

    % ------------------------------------------------------------------ %
    %  7. Pack raw result                                                  %
    % ------------------------------------------------------------------ %
    raw                  = struct();
    raw.B                = B_coef;          % [K*pmax x K] OLS coefficients
    raw.B_selection      = B;               % [K*pmax x K] binary mask from mBTS
    raw.Y0               = Y0;
    raw.residuals        = residuals;
    raw.fitted_values    = fitted;
    raw.T                = T;
    raw.K                = K;
    raw.T_eff            = T - pmax;
    raw.p  = pmax;   % framework uses pmax as the lag dimension of B
    raw.lambda           = NaN;             % mBTS is criterion-based, no lambda
    raw.toolbox_used     = 'none';

    % selected_terms: human-readable labels
    raw.selected_terms   = i_selected_term_labels(B_coef, K, pmax);

    % diagnostics
    raw.diagnostics = struct( ...
    'sparsity',         sum(B_coef(:) == 0) / numel(B_coef), ...
    'n_nonzero',        sum(B_coef(:) ~= 0), ...
    'mse_per_eq',       mean(residuals.^2, 1), ...
    'oos_msfe',         NaN, ...
    'maxorder_per_eq',  maxorders, ...
    'mse_per_eq_mbts',  mse_per_eq, ...
    'index_matrix',     indexM, ...
    'pmax',             pmax, ...
    'n_selected_total', sum(B(:)));

    raw.warnings = warn_msgs;

    raw.hyperparameters = struct( ...
        'pmax',    pmax, ...
        'method',  'mBTS', ...
        'criterion','BIC');

    % ------------------------------------------------------------------ %
    %  8. Normalize to canonical output                                   %
    % ------------------------------------------------------------------ %
    out = normalize_method_output(raw, METHOD, cfg, t_start);

    % ------------------------------------------------------------------ %
    %  9. mBTS-specific post-validation                                   %
    % ------------------------------------------------------------------ %
    out = i_mbts_postvalidate(out, K, pmax);
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function [Y0, X] = i_build_lagged(Y, T, K, pmax)
    % Builds lag-major regressor matrix consistent with mBTS indexV order.
    % X columns: [y1(t-1) y2(t-1)...yK(t-1) y1(t-2)...yK(t-pmax)]
    Y0 = Y(pmax+1:end, :);
    X  = zeros(T-pmax, K*pmax);
    for lag = 1:pmax
        cols       = (lag-1)*K + (1:K);
        X(:, cols) = Y(pmax+1-lag:end-lag, :);
    end
end

function labels = i_selected_term_labels(B, K, pmax)
    % Readable labels for non-zero entries of B selection mask
    [rows, cols] = find(B ~= 0);
    labels = cell(numel(rows), 1);
    for i = 1:numel(rows)
        lag_num    = ceil(rows(i) / K);
        var_num    = mod(rows(i)-1, K) + 1;
        eq_num     = cols(i);
        labels{i}  = sprintf('eq%d<-y%d_lag%d', eq_num, var_num, lag_num);
    end
end

function out = i_mbts_postvalidate(out, K, pmax)
    % mBTS-specific checks on top of generic validation
    if ~out.success
        return
    end

    % Warn if no terms selected anywhere
    if isfield(out, 'diagnostics') && ...
            isfield(out.diagnostics, 'n_selected_total') && ...
            out.diagnostics.n_selected_total == 0
        out.warnings{end+1} = ...
            'mBTS selected zero terms across all equations — check pmax or data.';
        out.status_code = 2;
    end

    % Warn if maxorder equals pmax in any equation (search may be truncated)
    if isfield(out.diagnostics, 'maxorder_per_eq')
        if any(out.diagnostics.maxorder_per_eq >= pmax)
            out.warnings{end+1} = sprintf( ...
                'maxorder reached pmax=%d in at least one equation — consider increasing pmax.', ...
                pmax);
            out.status_code = 2;
        end
    end
end

function cfg = i_fill_defaults(cfg)
    defaults = mbts_var_config();
    fields   = fieldnames(defaults);
    for i = 1:numel(fields)
        f = fields{i};
        if ~isfield(cfg, f)
            cfg.(f) = defaults.(f);
        end
    end
    % cfg.p and cfg.pmax sync
    if isempty(cfg.p)
        cfg.p = cfg.pmax;
    end
end

function out = i_early_failure(method, cfg, t_start, msg)
    raw         = struct();
    raw.success = false;
    raw.message = msg;
    out         = normalize_method_output(raw, method, cfg, t_start);
end