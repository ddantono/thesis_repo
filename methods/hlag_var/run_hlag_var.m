% methods/hlag_var/run_hlag_var.m

function out = run_hlag_var(Y, cfg)
% RUN_HLAG_VAR  HLAG-VAR via R/BigVAR (HLAGOO or HLAGC).
%
%   Hierarchical Lag (HLAG) structured penalty VAR, estimated via the
%   BigVAR R package. Two penalty structures are supported:
%
%     cfg.struct = 'HLAGOO'  — Own/Other hierarchical grouping
%     cfg.struct = 'HLAGC'   — Componentwise hierarchical grouping
%
%   Both share this wrapper. Dispatch keys: 'hlag_oo_var', 'hlag_c_var'.
%   The method field in cfg determines which is routed here; the struct
%   field tells the R script which BigVAR penalty to use.
%
%   Requires: R on system PATH, BigVAR R package.
%
%   Reference: Nicholson, Wilms, Bien & Matteson, JMLR 21(166), 2020.
%
%   INPUT
%     Y   : [T x K] double — stationary multivariate time series
%     cfg : struct from hlag_var_config()
%
%   OUTPUT
%     out : canonical struct (see create_canonical_output)

    t_start = tic;
    METHOD  = cfg.method;   % 'hlag_oo_var' or 'hlag_c_var'

    % ------------------------------------------------------------------ %
    %  1. Config                                                           %
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(cfg)
        cfg = hlag_var_config();
    end
    cfg = i_fill_defaults(cfg);

    % ------------------------------------------------------------------ %
    %  2. Validate struct field                                            %
    % ------------------------------------------------------------------ %
    valid_structs = {'HLAGOO', 'HLAGC'};
    if ~ismember(cfg.struct, valid_structs)
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('cfg.struct must be HLAGOO or HLAGC, got: %s', cfg.struct));
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

    % BigVAR needs enough obs for its rolling CV split
    % Minimum: T > 3*p + n_folds (heuristic guard)
    if T < 3*p + cfg.n_folds
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('T=%d too small for rolling CV with p=%d, n_folds=%d.', ...
                    T, p, cfg.n_folds));
        return
    end

    % ------------------------------------------------------------------ %
    %  4. Auto-compute T1, T2 if not supplied                             %
    % ------------------------------------------------------------------ %
    T1 = cfg.T1;
    T2 = cfg.T2;
    if isempty(T1) || isnan(T1), T1 = floor(T / 3);     end
    if isempty(T2) || isnan(T2), T2 = floor(2 * T / 3); end

    % ------------------------------------------------------------------ %
    %  5. Check R                                                          %
    % ------------------------------------------------------------------ %
    [r_ok, r_msg] = i_check_r(cfg.r_exe);
    if ~r_ok
        out = i_early_failure(METHOD, cfg, t_start, r_msg);
        return
    end

    % ------------------------------------------------------------------ %
    %  6. Locate R script                                                  %
    % ------------------------------------------------------------------ %
    r_script = i_resolve_r_script(cfg);
    if isempty(r_script)
        out = i_early_failure(METHOD, cfg, t_start, ...
            'hlag_var.R not found. Check cfg.r_script_path.');
        return
    end

    % ------------------------------------------------------------------ %
    %  7. Temp file paths                                                  %
    % ------------------------------------------------------------------ %
    run_id   = sprintf('%d', round(now * 1e6));
    in_file  = fullfile(cfg.tmp_dir, sprintf('hlag_input_%s.csv',  run_id));
    cfg_file = fullfile(cfg.tmp_dir, sprintf('hlag_config_%s.csv', run_id));
    out_file = fullfile(cfg.tmp_dir, sprintf('hlag_output_%s.csv', run_id));

    % ------------------------------------------------------------------ %
    %  8. Write input files                                                %
    % ------------------------------------------------------------------ %
    try
        write_r_input_file(Y, in_file);

        cfg_r               = struct();
        cfg_r.p             = p;
        cfg_r.K             = K;
        cfg_r.T             = T;
        cfg_r.struct        = cfg.struct;
        cfg_r.nlambda       = cfg.nlambda;
        cfg_r.n_folds       = cfg.n_folds;
        cfg_r.T1            = T1;
        cfg_r.T2            = T2;
        cfg_r.RVAR          = cfg.RVAR;
        cfg_r.Minnesota     = cfg.Minnesota;
        cfg_r.intercept     = cfg.intercept;
        cfg_r.rel_threshold = cfg.rel_threshold;   % ← προσθήκη

        write_r_config_file(cfg_r, cfg_file);

    catch ME
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Failed to write input files: %s', ME.message));
        return
    end

    % ------------------------------------------------------------------ %
    %  9. Call R                                                           %
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
    %  10. Read R output                                                   %
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
    %  11. Cleanup                                                         %
    % ------------------------------------------------------------------ %
    i_cleanup(cfg, in_file, cfg_file, out_file);

    % ------------------------------------------------------------------ %
    %  12. Normalize → canonical                                           %
    % ------------------------------------------------------------------ %
    out = normalize_r_output(r_result, METHOD, cfg, t_start, T, K, p, Y);
    
    if isfield(out.hyperparameters, 'alpha')
        out.hyperparameters = rmfield(out.hyperparameters, 'alpha');
    end
    if isfield(out.diagnostics, 'pilot')
        out.diagnostics = rmfield(out.diagnostics, 'pilot');
    end
    out.metadata.toolbox_used = 'R:BigVAR';

    % ------------------------------------------------------------------ %
    %  13. Method-specific post-validation                                 %
    % ------------------------------------------------------------------ %
    out = i_hlag_postvalidate(out, cfg);
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function cfg = i_fill_defaults(cfg)
    defaults = hlag_var_config();
    % Preserve caller's method field — do NOT overwrite with default
    caller_method = cfg.method;
    fields        = fieldnames(defaults);
    for i = 1:numel(fields)
        f = fields{i};
        if ~isfield(cfg, f)
            cfg.(f) = defaults.(f);
        end
    end
    cfg.method = caller_method;
    if isempty(cfg.p), cfg.p = 1; end
end

function [ok, msg] = i_check_r(r_exe)
    [status, ~] = system(sprintf('"%s" --version', r_exe));
    if status == 0
        ok  = true;
        msg = '';
    else
        ok  = false;
        msg = sprintf('R not found. Ensure "%s" is on system PATH.', r_exe);
    end
end

function r_script = i_resolve_r_script(cfg)
    if ~isempty(cfg.r_script_path) && exist(cfg.r_script_path, 'file')
        r_script = cfg.r_script_path;
        return
    end
    this_dir  = fileparts(mfilename('fullpath'));
    candidate = fullfile(this_dir, 'hlag_var.R');
    if exist(candidate, 'file')
        r_script = candidate;
    else
        r_script = '';
    end
end

function out = i_hlag_postvalidate(out, cfg)
    if ~out.success; return; end
    % Warn if all coefficients zero
    if ~isempty(out.coefficients) && all(out.coefficients(:) == 0)
        out.warnings{end+1} = ...
            'All coefficients zero — model may be over-regularised.';
        out.status_code = 2;
    end
    % Confirm struct stored correctly in hyperparameters
    if ~isfield(out.hyperparameters, 'struct') || ...
            ~strcmp(out.hyperparameters.struct, cfg.struct)
        out.warnings{end+1} = ...
            sprintf('hyperparameters.struct mismatch (expected %s).', cfg.struct);
        out.status_code = 2;
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