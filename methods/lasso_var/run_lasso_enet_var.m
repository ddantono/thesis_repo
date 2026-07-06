function out = run_lasso_enet_var(Y, cfg)
% RUN_LASSO_ENET_VAR  LASSO-VAR and Elastic Net-VAR wrapper.
%
%   cfg.alpha = 1.0        → LASSO
%   cfg.alpha ∈ (0, 1)     → Elastic Net
%   Estimates a sparse VAR(p) equation-by-equation using MATLAB's
%   built-in lasso(), which handles both LASSO (alpha=1) and
%   Elastic Net (0 < alpha < 1).
%
%   Requires: Statistics and Machine Learning Toolbox (lasso)
%
%   INPUT
%     Y   : [T x K] double — stationary multivariate time series
%     cfg : struct from lasso_var_config() or user-supplied
%
%   OUTPUT
%     out : canonical struct (see create_canonical_output)

    t_start = tic;
    METHOD  = 'lasso_enet_var';

    % ------------------------------------------------------------------ %
    %  1. Config fallback                                                  %
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(cfg)
        cfg = lasso_var_config();
    end
    cfg = i_fill_defaults(cfg);

    % ------------------------------------------------------------------ %
    %  2. Toolbox check                                                    %
    % ------------------------------------------------------------------ %
    if ~license('test', 'statistics_toolbox')
        out = i_early_failure(METHOD, cfg, t_start, ...
            'Statistics and Machine Learning Toolbox not available.');
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
    if isfield(cfg, 'pmax')
        p = cfg.pmax;
    elseif isfield(cfg, 'p')
        p = cfg.p;
    else
        error('run_lasso_enet_var: no lag order found in cfg');
    end

    if T <= p + 1
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Insufficient observations T=%d for p=%d.', T, p));
        return
    end

    % ------------------------------------------------------------------ %
    %  4. Build lagged regressor matrix X and response Y0                 %
    % ------------------------------------------------------------------ %
    [Y0, X] = i_build_lagged(Y, T, K, p);   % Y0: [(T-p) x K], X: [(T-p) x K*p]
    T_eff   = T - p;

    % ------------------------------------------------------------------ %
    %  5. Equation-by-equation LASSO / Elastic Net                        %
    % ------------------------------------------------------------------ %
    B          = zeros(K*p, K);
    lambda_sel = zeros(1, K);       % selected lambda per equation
    fit_info   = cell(1, K);        % lasso() FitInfo structs
    warn_msgs  = {};

    lasso_opts = i_build_lasso_opts(cfg);

    for k = 1:K
        y_k = Y0(:, k);

        try
            if cfg.cv
                % CV-based lambda selection via cvpartition
                [B_path, FitInfo] = lasso(X, y_k, lasso_opts{:});
                idx_sel           = i_select_lambda(FitInfo, cfg.cv_criterion);
                B(:, k)           = B_path(:, idx_sel);
                lambda_sel(k)     = FitInfo.Lambda(idx_sel);
            else
                % Fixed lambda supplied by user
                [B_path, FitInfo] = lasso(X, y_k, lasso_opts{:}, ...
                                          'Lambda', cfg.lambda);
                B(:, k)           = B_path(:, 1);
                lambda_sel(k)     = cfg.lambda;
            end
            fit_info{k} = FitInfo;

        catch ME
            warn_msgs{end+1} = sprintf('Equation %d failed: %s', k, ME.message); %#ok<AGROW>
            % B(:,k) stays zero — flag but continue
        end
    end

    % ------------------------------------------------------------------ %
    %  6. Derived quantities                                               %
    % ------------------------------------------------------------------ %
    fitted    = X * B;                  % [(T-p) x K]
    residuals = Y0 - fitted;            % [(T-p) x K]

    % ------------------------------------------------------------------ %
    %  7. Pack raw_result for normalize_method_output                     %
    % ------------------------------------------------------------------ %
    raw               = struct();
    raw.B             = B;
    raw.Y0            = Y0;
    raw.residuals     = residuals;
    raw.fitted_values = fitted;
    raw.T             = T;
    raw.K             = K;
    raw.T_eff         = T_eff;
    raw.p             = p;
    raw.lambda        = lambda_sel;   % [1 x K] one per equation
    raw.toolbox_used  = 'statistics_toolbox:lasso';

    % selected_terms: cell of labels for non-zero coefficient rows
    raw.selected_terms = i_selected_term_labels(B, K, p);

    % method-specific diagnostics
    raw.diagnostics = i_build_diagnostics(B, residuals, lambda_sel, ...
                                          fit_info, cfg, K, p);

    % warnings from failed equations
    raw.warnings = warn_msgs;

    % hyperparameters actually used
    raw.hyperparameters = struct( ...
        'alpha',       cfg.alpha, ...
        'lambda',      lambda_sel, ...
        'cv',          cfg.cv, ...
        'n_folds',     cfg.n_folds, ...
        'cv_criterion',cfg.cv_criterion, ...
        'standardize', cfg.standardize);

    % ------------------------------------------------------------------ %
    %  8. Normalize to canonical output                                    %
    % ------------------------------------------------------------------ %
    out = normalize_method_output(raw, METHOD, cfg, t_start);

    % ------------------------------------------------------------------ %
    %  9. Method-specific post-validation                                  %
    % ------------------------------------------------------------------ %
    out = i_lasso_postvalidate(out, cfg);
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function [Y0, X] = i_build_lagged(Y, T, K, p)
    Y0 = Y(p+1:end, :);
    X  = zeros(T-p, K*p);
    for lag = 1:p
        cols       = (lag-1)*K + (1:K);
        X(:, cols) = Y(p+1-lag:end-lag, :);
    end
end

function opts = i_build_lasso_opts(cfg)
    % Build name-value pairs for lasso() that are always passed
    opts = { ...
        'Alpha',       cfg.alpha, ...
        'Standardize', cfg.standardize, ...
        'MaxIter',     cfg.max_iter, ...
        'RelTol',      cfg.rel_tol};

    if cfg.cv
        opts = [opts, {'CV', cfg.n_folds}];
    end
end

function idx = i_select_lambda(FitInfo, criterion)
    % Select lambda index from CV FitInfo
    switch lower(criterion)
        case 'min+1se'
            idx = FitInfo.Index1SE;   % lasso() built-in 1-SE rule
        otherwise
            idx = FitInfo.IndexMinMSE;
    end
    if isempty(idx) || isnan(idx)
        idx = 1;   % fallback: most regularised
    end
end

function labels = i_selected_term_labels(B, K, p)
    % Generate readable labels for non-zero rows of B
    nz     = any(B ~= 0, 2);
    idx    = find(nz);
    labels = cell(numel(idx), 1);
    for i = 1:numel(idx)
        lag_num  = ceil(idx(i) / K);
        var_num  = mod(idx(i)-1, K) + 1;
        labels{i} = sprintf('y%d_lag%d', var_num, lag_num);
    end
end

function diag = i_build_diagnostics(B, residuals, lambda_sel, ...
                                     fit_info, cfg, K, p)
    diag.sparsity     = sum(B(:) == 0) / numel(B);
    diag.n_nonzero    = sum(B(:) ~= 0);
    diag.mse_per_eq   = mean(residuals.^2, 1);
    diag.oos_msfe     = NaN;
    diag.lambda_per_eq= lambda_sel;
    diag.alpha        = cfg.alpha;
    diag.method_type  = i_method_type(cfg.alpha);
    diag.sigma        = (residuals' * residuals) / size(residuals, 1);

    % Extract CV MSE path if available
    if cfg.cv && ~isempty(fit_info{1}) && ...
            isfield(fit_info{1}, 'MSE')
        diag.cv_mse_path = cellfun(@(f) f.MSE, fit_info, ...
                                   'UniformOutput', false);
    else
        diag.cv_mse_path = {};
    end
end

function s = i_method_type(alpha)
    if alpha == 1
        s = 'LASSO';
    elseif alpha > 0 && alpha < 1
        s = 'ElasticNet';
    else
        s = 'Ridge';   % alpha=0 edge case
    end
end

function out = i_lasso_postvalidate(out, cfg)
    % LASSO-specific checks on top of generic validation
    if out.success
        % Warn if all coefficients are zero (over-regularised)
        if ~isempty(out.coefficients) && all(out.coefficients(:) == 0)
            out.warnings{end+1} = ...
                'All coefficients are zero — model may be over-regularised.';
            out.status_code = 2;
        end

        % Warn if alpha is outside (0,1]
        if cfg.alpha <= 0 || cfg.alpha > 1
            out.warnings{end+1} = ...
                sprintf('alpha=%.4f is outside (0,1]. Check config.', cfg.alpha);
            out.status_code = 2;
        end
    end
end

function cfg = i_fill_defaults(cfg)
    defaults = lasso_var_config();
    fields   = fieldnames(defaults);
    for i = 1:numel(fields)
        f = fields{i};
        if ~isfield(cfg, f)
            cfg.(f) = defaults.(f);
        end
    end
end

function out = i_early_failure(method, cfg, t_start, msg)
    raw         = struct();
    raw.success = false;
    raw.message = msg;
    out         = normalize_method_output(raw, method, cfg, t_start);
end