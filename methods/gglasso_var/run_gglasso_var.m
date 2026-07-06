% methods/gglasso_var/run_gglasso_var.m

function out = run_gglasso_var(Y, cfg)
% RUN_GGLASSO_VAR  Group LASSO VAR via R/gglasso.
%
%   Fits a sparse VAR(p) equation-by-equation using group lasso penalty.
%   Groups are formed lag-wise (default): all K predictor coefficients at
%   a given lag for one equation constitute one group.  This promotes
%   entire-lag sparsity (either all variables at lag l enter or none do).
%
%   Alternatively, groups can be formed variable-wise (cfg.group_by =
%   'variable'): all p lags of one predictor form a group, promoting
%   predictor-level exclusion across lags.
%
%   CV (k-fold, criterion lambda.min) selects lambda per equation.
%
%   Requires: R on PATH, gglasso R package.
%   Reference: Yang & Zou (2015), Statistics and Computing 25(6).
%
%   INPUT
%     Y   : [T x K] double
%     cfg : struct from gglasso_var_config()
%
%   OUTPUT
%     out : canonical struct (see create_canonical_output)

    t_start = tic;
    METHOD  = 'gglasso_var';

    % ------------------------------------------------------------------ %
    %  1. Config                                                           %
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(cfg)
        cfg = gglasso_var_config();
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

    if T <= p + 1
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Insufficient observations T=%d for p=%d.', T, p));
        return
    end
    if K < 2
        out = i_early_failure(METHOD, cfg, t_start, ...
            'Group LASSO VAR requires K >= 2.');
        return
    end

    % ------------------------------------------------------------------ %
    %  3. Check R availability                                             %
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
            'gglasso_var.R not found. Check cfg.r_script_path.');
        return
    end

    % ------------------------------------------------------------------ %
    %  5. Prepare temp file paths                                          %
    % ------------------------------------------------------------------ %
    run_id   = sprintf('%d', round(now * 1e6));
    in_file  = fullfile(cfg.tmp_dir, sprintf('gglasso_input_%s.csv',  run_id));
    out_file = fullfile(cfg.tmp_dir, sprintf('gglasso_output_%s.csv', run_id));
    cfg_file = fullfile(cfg.tmp_dir, sprintf('gglasso_config_%s.csv', run_id));

    % ------------------------------------------------------------------ %
    %  6. Write input files                                                %
    % ------------------------------------------------------------------ %
    try
        write_r_input_file(Y, in_file);

        cfg_r              = struct();
        cfg_r.p            = p;
        cfg_r.K            = K;
        cfg_r.nlambda      = cfg.nlambda;
        cfg_r.loss         = cfg.loss;
        cfg_r.cv_criterion = cfg.cv_criterion;
        cfg_r.n_folds      = cfg.n_folds;
        cfg_r.eps          = cfg.eps;
        cfg_r.maxit        = cfg.maxit;
        cfg_r.intercept    = cfg.intercept;
        cfg_r.group_by     = cfg.group_by;
        % lambda: -1 sentinel means auto-select in R
        if isempty(cfg.lambda)
            cfg_r.lambda   = -1;
        else
            cfg_r.lambda   = cfg.lambda;
        end

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
    %  9. Cleanup temp files                                               %
    % ------------------------------------------------------------------ %
    i_cleanup(cfg, in_file, cfg_file, out_file);

    % ------------------------------------------------------------------ %
    %  10. Normalize to canonical output                                  %
    % ------------------------------------------------------------------ %
    out = normalize_r_output(r_result, METHOD, cfg, t_start, T, K, p, Y);

    % ------------------------------------------------------------------ %
    %  11. Fix generic fields that normalize_r_output sets incorrectly    %
    %      for gglasso (no glmnet, no pilot, no alpha in the glmnet sense)%
    % ------------------------------------------------------------------ %

    % alpha in gglasso context is not the elastic-net mixing parameter;
    % remove it from hyperparameters to avoid confusion
    if isfield(out.hyperparameters, 'alpha')
        out.hyperparameters = rmfield(out.hyperparameters, 'alpha');
    end

    % pilot is a glmnet/adaptive-lasso concept — not applicable here
    if isfield(out.diagnostics, 'pilot')
        out.diagnostics = rmfield(out.diagnostics, 'pilot');
    end

    % toolbox label
    out.metadata.toolbox_used = 'R:gglasso';

    % group_by scheme into hyperparameters
    out.hyperparameters.group_by = cfg.group_by;
    out.hyperparameters.loss     = cfg.loss;

    % oos_msfe: R writes literal "NaN" string in some paths → scalar NaN
    if ischar(out.diagnostics.oos_msfe) || isstring(out.diagnostics.oos_msfe)
        out.diagnostics.oos_msfe = NaN;
    end

    % n_groups from diagnostics (if R wrote it) — move to hyperparameters
    if isfield(out.diagnostics, 'n_groups')
        out.hyperparameters.n_groups = out.diagnostics.n_groups;
        out.diagnostics = rmfield(out.diagnostics, 'n_groups');
    end
    
    % Verify centering was effective: check R-side intercept values
    % If centering worked, intercept_vals should be ~0
    if ~cfg.intercept && isfield(out.hyperparameters, 'intercept')
        intc = out.hyperparameters.intercept;
        if isnumeric(intc) && any(abs(intc) > 1e-4)
            out.warnings{end+1} = sprintf( ...
                ['gglasso: intercept=false with R-side centering, but ' ...
                 'fitted intercept max(|a|)=%.4g — centering may have failed.'], ...
                max(abs(intc)));
            if out.status_code == 0
                out.status_code = 2;
            end
        end
    end

    % ------------------------------------------------------------------ %
    %  12. Method-specific post-validation                                %
    % ------------------------------------------------------------------ %

    out = i_gglasso_postvalidate(out, cfg);
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function cfg = i_fill_defaults(cfg)
    defaults = gglasso_var_config();
    fields   = fieldnames(defaults);
    for i = 1:numel(fields)
        f = fields{i};
        if ~isfield(cfg, f)
            cfg.(f) = defaults.(f);
        end
    end
    if isempty(cfg.p) || cfg.p < 1
        cfg.p = 1;
    end
end

function [ok, msg] = i_check_r(r_exe)
    [status, ~] = system(sprintf('"%s" --version', r_exe));
    if status == 0
        ok  = true;
        msg = '';
    else
        ok  = false;
        msg = sprintf('R not found at "%s". Ensure R is on system PATH.', r_exe);
    end
end

function r_script = i_resolve_r_script(cfg)
    % 1. Explicit path
    if ~isempty(cfg.r_script_path) && exist(cfg.r_script_path, 'file')
        r_script = cfg.r_script_path;
        return
    end
    % 2. Sibling to this wrapper (requires methods/gglasso_var/ on MATLAB path)
    wrapper_path = which('run_gglasso_var');
    if ~isempty(wrapper_path)
        candidate = fullfile(fileparts(wrapper_path), 'gglasso_var.R');
        if exist(candidate, 'file')
            r_script = candidate;
            return
        end
    end
    r_script = '';
end

function out = i_gglasso_postvalidate(out, cfg)
    if ~out.success; return; end

    % Warn if all coefficients zero — model over-regularised
    if ~isempty(out.coefficients) && all(out.coefficients(:) == 0)
        out.warnings{end+1} = ...
            'All group-lasso coefficients zero — model may be over-regularised.';
        out.status_code = 2;
    end

    % Warn if CV selected maximum lambda (edge of path)
    if isfield(out.hyperparameters, 'lambda') && ...
            isnumeric(out.hyperparameters.lambda)
        lam = out.hyperparameters.lambda;
        if numel(lam) == 1 && ~isnan(lam) && lam > 0 && ...
                isfield(out.diagnostics, 'lambda_max') && ...
                lam >= out.diagnostics.lambda_max * 0.99
            out.warnings{end+1} = ...
                'CV selected lambda is at the top of the path — all groups may be zero.';
            out.status_code = 2;
        end
    end

    % Granger causality adjacency from group-lasso B
    K = out.metadata.K;
    p = out.lag_order;
    if ~isempty(out.coefficients) && ~isnan(K) && ~isnan(p)
        B   = out.coefficients;   % [K*p x K]
        adj = zeros(K, K);
        for cause = 1:K
            rows = ((0:p-1)*K) + cause;
            rows = rows(rows <= size(B,1));
            for eq = 1:K
                if cause ~= eq
                    adj(eq, cause) = any(B(rows, eq) ~= 0);
                end
            end
        end
        out.diagnostics.gc_matrix = adj;   % gc(i,j)=1 → yj Granger-causes yi
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