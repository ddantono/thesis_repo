% methods/msvar/run_msvar.m

function out = run_msvar(Y, cfg)
% RUN_MSVAR  Two-stage modified sparse VAR (msVAR) via R/TSGlasso.
%
%   Stage 1: TSGlasso identifies non-zero pairs of the inverse spectral
%            density matrix; BIC selects the lag order p.
%   Stage 2: FDR (Benjamini-Hochberg) refines non-zero AR coefficients.
%
%   The R script (msvar.R) sources the two supplementary files:
%     msVAR_function.r  (mmc2) — core msVAR machinery
%     tsGLASSO.R        (mmc3) — TSGlasso ADMM solver
%   Both must reside in the same directory as msvar.R.
%
%   Requires: R on PATH, LongMemoryTS R package.
%   Reference: Dallakyan, Kim & Pourahmadi (2022), CSDA 176, 107557.
%
%   INPUT
%     Y   : [T x K] double — stationary multivariate time series
%     cfg : struct from msvar_config()
%
%   OUTPUT
%     out : canonical struct (see create_canonical_output)

    t_start = tic;
    METHOD  = 'msvar';

    % ------------------------------------------------------------------ %
    %  1. Config                                                           %
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(cfg)
        cfg = msvar_config();
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
    p_max  = max(cfg.p_seq);

    if T <= p_max + K + 1
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Insufficient observations T=%d for p_max=%d, K=%d.', ...
                    T, p_max, K));
        return
    end
    if K < 3
        out = i_early_failure(METHOD, cfg, t_start, ...
            'msVAR requires K >= 3 (TSGlasso constraint).');
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
    %  4. Locate R script and its dependencies                             %
    % ------------------------------------------------------------------ %
    [r_script, dep_msg] = i_resolve_r_script(cfg);
    if isempty(r_script) || ~isempty(dep_msg)
        out = i_early_failure(METHOD, cfg, t_start, dep_msg);
        return
    end

    % ------------------------------------------------------------------ %
    %  5. Prepare temp file paths                                          %
    % ------------------------------------------------------------------ %
    run_id   = sprintf('%d', round(now * 1e6));
    in_file  = fullfile(cfg.tmp_dir, sprintf('msvar_input_%s.csv',  run_id));
    out_file = fullfile(cfg.tmp_dir, sprintf('msvar_output_%s.csv', run_id));
    cfg_file = fullfile(cfg.tmp_dir, sprintf('msvar_config_%s.csv', run_id));

    % ------------------------------------------------------------------ %
    %  6. Write input files                                                %
    % ------------------------------------------------------------------ %
    try
        write_r_input_file(Y, in_file);

        cfg_r                      = struct();
        cfg_r.K                    = K;
        % p_seq written as a single string because write_r_config_file
        % cannot serialise vectors — R parses it with scan()
        cfg_r.p_seq_str            = strjoin(arrayfun(@num2str, cfg.p_seq, ...
                                       'UniformOutput', false), ' ');
        if isempty(cfg.lambda)
            cfg_r.lambda           = -1;    % sentinel: auto-select in R
        else
            cfg_r.lambda           = cfg.lambda;
        end
        cfg_r.half_window_length   = cfg.half_window_length;
        cfg_r.rho                  = cfg.rho;
        cfg_r.alpha                = cfg.alpha;
        cfg_r.thresh               = cfg.thresh;
        cfg_r.rho_flex             = cfg.rho_flex;
        cfg_r.ADMM_iter            = cfg.ADMM_iter;
        cfg_r.stage1_info          = cfg.stage1_info;
        cfg_r.standardize          = cfg.standardize;
        cfg_r.lambda_gam           = cfg.lambda_gam;
        cfg_r.lambda_trim_max      = cfg.lambda_trim_max;
        cfg_r.lambda_trim_min      = cfg.lambda_trim_min;
        cfg_r.lambda_criteria      = cfg.lambda_criteria;
        cfg_r.fdr_q                = cfg.fdr_q;
        cfg_r.ite_max              = cfg.ite_max;
        cfg_r.reltol               = cfg.reltol;

        write_r_config_file(cfg_r, cfg_file);

    catch ME
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Failed to write input files: %s', ME.message));
        return
    end

    % ------------------------------------------------------------------ %
    %  7. Call R                                                           %
    %  Pass script_dir explicitly so R can source() its dependencies      %
    %  regardless of working directory.                                    %
    % ------------------------------------------------------------------ %
    script_dir = fileparts(r_script);   % directory containing msvar.R
    [r_ok, r_log] = call_r_script(cfg.r_exe, r_script, ...
                                   in_file, cfg_file, script_dir, out_file);
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
    %  10. Determine p values from R output                               %
    %                                                                     %
    %  p_selected : BIC-selected lag (controls gc_matrix, postvalidate)  %
    %  p_ub       : max(p_seq) — B matrix always has K*p_ub rows         %
    %  normalize_r_output must receive p_ub so X*B dimensions match.     %
    % ------------------------------------------------------------------ %
    p_selected = i_extract_p(r_result, cfg.p_seq);
    p_ub       = max(cfg.p_seq);

    % cfg.p = p_ub so framework dimension checks pass
    cfg.p = p_ub;

    % ------------------------------------------------------------------ %
    %  11. Normalize to canonical output                                   %
    % ------------------------------------------------------------------ %
    out = normalize_r_output(r_result, METHOD, cfg, t_start, T, K, p_ub, Y);

    % Stamp the BIC-selected lag in lag_order and hyperparameters
    out.lag_order                        = p_selected;
    % Truncate B to K*p_selected rows (extra rows are zero by construction)
    if p_selected < p_ub && ~isempty(out.coefficients)
        out.coefficients = out.coefficients(1:K*p_selected, :);
        % Update residuals/fitted from truncated B
        Y0 = Y(p_selected+1:end, :);
        X  = zeros(T-p_selected, K*p_selected);
        for lag = 1:p_selected
            cols = (lag-1)*K + (1:K);
            X(:, cols) = Y(p_selected+1-lag:end-lag, :);
        end
        out.residuals    = Y0 - X * out.coefficients;
        out.fitted_values = X * out.coefficients;
        out.metadata.T_eff = T - p_selected;
    end
    out.hyperparameters.p_selected       = p_selected;
    out.hyperparameters.p_ub             = p_ub;
    out.metadata.toolbox_used = 'R:custom';
    
    % ------------------------------------------------------------------ %
    %  Post-processing: fix fields that normalize_r_output sets generically%
    % ------------------------------------------------------------------ %

    % 1. oos_msfe: R writes literal "NaN" string → convert to scalar NaN
    if ischar(out.diagnostics.oos_msfe) || isstring(out.diagnostics.oos_msfe)
        out.diagnostics.oos_msfe = NaN;
    end

    % 2. alpha: normalize_r_output pulls cfg.alpha = ADMM relaxation (1.5)
    %    msVAR does not use glmnet alpha — remove to avoid confusion
    if isfield(out.hyperparameters, 'alpha')
        out.hyperparameters = rmfield(out.hyperparameters, 'alpha');
    end

    % 3. mse_per_eq in hyperparameters: duplicate of diagnostics — remove
    if isfield(out.hyperparameters, 'mse_per_eq')
        out.hyperparameters = rmfield(out.hyperparameters, 'mse_per_eq');
    end

    % 4. cv_criterion / n_folds / standardize: not applicable to msVAR
    for f = {'cv_criterion', 'n_folds', 'standardize'}
        if isfield(out.hyperparameters, f{1})
            out.hyperparameters = rmfield(out.hyperparameters, f{1});
        end
    end

    % 5b. pilot in diagnostics: artifact of normalize_r_output — remove
    if isfield(out.diagnostics, 'pilot')
        out.diagnostics = rmfield(out.diagnostics, 'pilot');
    end

    % 6. Intercept diagnostics — pulled from R meta
    %    intercept_max   : max|intercept_k| — if large, data is not centered
    %    intercept_mean  : mean|intercept_k|
    %    intercept_centered : 1 if all intercepts < 1e-6 (effectively zero)
    %    These inform whether mse_per_eq is affected by omitted intercept.
    for f = {'intercept_max', 'intercept_mean', 'intercept_centered'}
        fname = f{1};
        if isfield(out.hyperparameters, fname)
            out.diagnostics.(fname) = out.hyperparameters.(fname);
            out.hyperparameters     = rmfield(out.hyperparameters, fname);
        end
    end

    % Intercept warning: trigger only when intercept materially affects MSE.
    % Criterion: max delta_mse% > 5% across equations.
    % (intercept/residual_std is not reliable — validated against manual check)
    if isfield(out.diagnostics, 'intercept_max') && ...
            isnumeric(out.diagnostics.intercept_max) && ...
            ~isnan(out.diagnostics.intercept_max) && ...
            out.diagnostics.intercept_max > 0 && ...
            ~isempty(out.residuals)
        mse_no_int  = mean(out.residuals.^2, 1);
        mse_with_in = out.diagnostics.mse_per_eq;
        if isnumeric(mse_with_in) && numel(mse_with_in) == size(out.residuals,2)
            delta_pct = 100 * abs(mse_with_in - mse_no_int) ./ (mse_no_int + eps);
            max_delta = max(delta_pct);
            % Store delta for Section 11 of eval scripts
            out.diagnostics.intercept_delta_mse_pct = max_delta;
            if max_delta > 5
                out.warnings{end+1} = sprintf( ...
                    ['msVAR: intercept causes %.1f%% MSE difference. ' ...
                     'Center data (Y = Y - mean(Y)) for fair comparison.'], ...
                    max_delta);
                if out.status_code == 0
                    out.status_code = 2;
                end
            end
        end
    end

    % 5. lambda_per_eq in diagnostics: msVAR has one global lambda
    %    move it to hyperparameters.lambda if not already there
    if isfield(out.diagnostics, 'lambda_per_eq')
        if ~isfield(out.hyperparameters, 'lambda') || ...
                isempty(out.hyperparameters.lambda)
            out.hyperparameters.lambda = out.diagnostics.lambda_per_eq;
        end
        out.diagnostics = rmfield(out.diagnostics, 'lambda_per_eq');
    end

    % ------------------------------------------------------------------ %
    %  12. Method-specific post-validation                                 %
    % ------------------------------------------------------------------ %
    out = i_msvar_postvalidate(out, K, p_selected, cfg);
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function cfg = i_fill_defaults(cfg)
    defaults = msvar_config();
    fields   = fieldnames(defaults);
    for i = 1:numel(fields)
        f = fields{i};
        if ~isfield(cfg, f) || isempty(cfg.(f))
            cfg.(f) = defaults.(f);
        end
    end
    % p_seq must be a non-empty sorted integer vector
    if ~isfield(cfg, 'p_seq') || isempty(cfg.p_seq)
        cfg.p_seq = [1 2 3];
    end
    cfg.p_seq = sort(unique(round(cfg.p_seq)));
    % cfg.p: use max(p_seq) as placeholder until R returns the BIC winner
    if isempty(cfg.p)
        cfg.p = max(cfg.p_seq);
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

function [r_script, msg] = i_resolve_r_script(cfg)
% Resolves the R script path.
% Priority:
%   1. cfg.r_script_path (explicit, always preferred)
%   2. which('run_msvar') → infer sibling msvar.R  (requires MATLAB path)
%   3. Fail with a clear message.

    msg = '';

    % --- 1. Explicit user-supplied path ----------------------------------
    if isfield(cfg, 'r_script_path') && ~isempty(cfg.r_script_path)
        if exist(cfg.r_script_path, 'file')
            r_script = cfg.r_script_path;
            msg      = i_check_deps(fileparts(r_script));
            return
        else
            r_script = '';
            msg = sprintf('cfg.r_script_path not found: %s', cfg.r_script_path);
            return
        end
    end

    % --- 2. Infer from which('run_msvar') --------------------------------
    wrapper_path = which('run_msvar');
    if ~isempty(wrapper_path)
        this_dir  = fileparts(wrapper_path);
        candidate = fullfile(this_dir, 'msvar.R');
        if exist(candidate, 'file')
            r_script = candidate;
            msg      = i_check_deps(this_dir);
            return
        end
    end

    % --- 3. Fail ---------------------------------------------------------
    r_script = '';
    msg = ['msvar.R not found. Set cfg.r_script_path to the full path of ' ...
           'msvar.R, or add methods/msvar/ to the MATLAB path.'];
end

function msg = i_check_deps(dir_path)
% Returns non-empty error message if companion R scripts are missing.
    msg  = '';
    deps = {'msVAR_function.r', 'tsGLASSO.r'};
    for i = 1:numel(deps)
        if ~exist(fullfile(dir_path, deps{i}), 'file')
            msg = sprintf(['Dependency "%s" not found in %s. ' ...
                'Place msVAR_function.r (mmc2) and tsGLASSO.R (mmc3) ' ...
                'in the same folder as msvar.R.'], deps{i}, dir_path);
            return
        end
    end
end

function p = i_extract_p(r_result, p_seq)
    % R writes the BIC-selected lag order into meta.p_selected
    p = NaN;
    if isfield(r_result, 'meta') && isfield(r_result.meta, 'p_selected')
        val = r_result.meta.p_selected;
        if isnumeric(val) && isscalar(val) && val >= 1
            p = round(val);
            return
        end
    end
    % Fallback: infer from B dimensions
    if isfield(r_result, 'B') && ~isempty(r_result.B)
        K_p = size(r_result.B, 1);
        K   = size(r_result.B, 2);
        if K > 0 && mod(K_p, K) == 0
            p = K_p / K;
            return
        end
    end
    % Last resort: use max(p_seq)
    p = max(p_seq);
end

function out = i_msvar_postvalidate(out, K, p, cfg)
    if ~out.success; return; end

    % Warn if all coefficients are zero
    if ~isempty(out.coefficients) && all(out.coefficients(:) == 0)
        out.warnings{end+1} = ...
            'All Stage-2 coefficients zero — model may be over-regularised.';
        out.status_code = 2;
    end

    % Warn if selected p is at the boundary of the search space
    if ~isnan(p) && p == max(cfg.p_seq)
        out.warnings{end+1} = sprintf( ...
            'BIC selected p=%d (maximum of search space). Consider expanding p_seq.', p);
        out.status_code = 2;
    end

    % Granger causality adjacency matrix derived from Stage-2 B
    if ~isempty(out.coefficients) && ~isnan(p) && ~isnan(K)
        B   = out.coefficients;   % [K*p x K]
        adj = zeros(K, K);
        for cause = 1:K
            rows = ((0:p-1) * K) + cause;   % rows of B for this cause variable
            rows = rows(rows <= size(B,1));
            for eq = 1:K
                if cause ~= eq
                    adj(eq, cause) = any(B(rows, eq) ~= 0);
                end
            end
        end
        out.diagnostics.gc_matrix = adj;   % [K x K], gc(i,j)=1 → yj→yi
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