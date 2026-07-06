% methods/ncpen_var/run_ncpen_var.m

function out = run_ncpen_var(Y, cfg)
% RUN_NCPEN_VAR  Nonconcave Penalized VAR (SCAD/MCP) via R/bigVAR.
%
%   Delegates estimation to R using the bigVAR package.
%   Supports SCAD and MCP penalties via cfg.penalty.
%
%   Requires: R on system PATH, bigVAR R package,
%             core/r_bridge/ functions on MATLAB path.
%
%   Reference: Nicholson, Matteson & Bien (2017), bigVAR R package.
%              Davis, Zang & Zheng (2016), JCGS.
%
%   INPUT
%     Y   : [T x K] double — stationary multivariate time series
%     cfg : struct from ncpen_var_config()
%
%   OUTPUT
%     out : canonical struct (see create_canonical_output)

    t_start = tic;
    METHOD  = 'ncpen_var';

    % ------------------------------------------------------------------ %
    %  1. Config                                                           %
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(cfg)
        cfg = ncpen_var_config();
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

    % Validate penalty
    if ~ismember(upper(cfg.penalty), {'SCAD', 'MCP'})
        out = i_early_failure(METHOD, cfg, t_start, ...
            sprintf('Invalid penalty "%s". Use "SCAD" or "MCP".', cfg.penalty));
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
            'ncpen_var.R not found. Check cfg.r_script_path.');
        return
    end

    % ------------------------------------------------------------------ %
    %  5. Prepare temp file paths                                          %
    % ------------------------------------------------------------------ %
    run_id   = sprintf('%d', round(now * 1e6));
    in_file  = fullfile(cfg.tmp_dir, sprintf('ncpen_input_%s.csv',  run_id));
    out_file = fullfile(cfg.tmp_dir, sprintf('ncpen_output_%s.csv', run_id));
    cfg_file = fullfile(cfg.tmp_dir, sprintf('ncpen_config_%s.csv', run_id));

    % ------------------------------------------------------------------ %
    %  6. Write input files                                                %
    % ------------------------------------------------------------------ %
    try
        write_r_input_file(Y, in_file);

        cfg_r             = struct();
        cfg_r.p           = p;
        cfg_r.K           = K;
        cfg_r.penalty     = upper(cfg.penalty);
        cfg_r.n_folds     = cfg.n_folds;
        cfg_r.nlambda     = cfg.nlambda;
        cfg_r.Minnesota   = cfg.Minnesota;
        cfg_r.verbose     = cfg.verbose;
        if ~isempty(cfg.lambda)
            cfg_r.lambda  = cfg.lambda;
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
    %  9. Cleanup                                                          %
    % ------------------------------------------------------------------ %
    i_cleanup(cfg, in_file, cfg_file, out_file);

    % ------------------------------------------------------------------ %
    %  10. Normalize to canonical output                                   %
    % ------------------------------------------------------------------ %
    out = normalize_r_output(r_result, METHOD, cfg, t_start, T, K, p, Y);

    % Fix generic fields from normalize_r_output (not glmnet-based)
    if isfield(out.hyperparameters, 'alpha')
        out.hyperparameters = rmfield(out.hyperparameters, 'alpha');
    end
    if isfield(out.diagnostics, 'pilot')
        out.diagnostics = rmfield(out.diagnostics, 'pilot');
    end
    out.metadata.toolbox_used = 'R:BigVAR';

    % ------------------------------------------------------------------ %
    %  11. Method-specific post-validation                                 %
    % ------------------------------------------------------------------ %
    out = i_ncpen_postvalidate(out, cfg);
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function cfg = i_fill_defaults(cfg)
    defaults = ncpen_var_config();
    fields   = fieldnames(defaults);
    for i = 1:numel(fields)
        f = fields{i};
        if ~isfield(cfg, f)
            cfg.(f) = defaults.(f);
        end
    end
    if isempty(cfg.p)
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
        msg = sprintf('R not found. Ensure "%s" is on system PATH.', r_exe);
    end
end

function r_script = i_resolve_r_script(cfg)
    if ~isempty(cfg.r_script_path) && exist(cfg.r_script_path, 'file')
        r_script = cfg.r_script_path;
        return
    end
    this_dir  = fileparts(mfilename('fullpath'));
    candidate = fullfile(this_dir, 'ncpen_var.R');
    if exist(candidate, 'file')
        r_script = candidate;
    else
        r_script = '';
    end
end

function out = i_ncpen_postvalidate(out, cfg)
    if ~out.success; return; end

    % Warn if all coefficients zero
    if ~isempty(out.coefficients) && all(out.coefficients(:) == 0)
        out.warnings{end+1} = ...
            'All coefficients zero — model may be over-regularised.';
        out.status_code = 2;
    end

    % Warn if penalty not recognised (defensive)
    if ~ismember(upper(cfg.penalty), {'SCAD', 'MCP'})
        out.warnings{end+1} = sprintf('Unrecognised penalty: %s.', cfg.penalty);
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