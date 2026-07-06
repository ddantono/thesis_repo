function out = normalize_method_output(raw, method_name, cfg, t_start)
% NORMALIZE_METHOD_OUTPUT  Converts any method's raw output to canonical form.
%
%   out = normalize_method_output(raw, method_name, cfg, t_start)
%
%   raw         : whatever the method returned (struct expected; scalar/[]
%                 accepted as failure sentinel)
%   method_name : char — must match dispatch key
%   cfg         : config struct used for the run
%   t_start     : tic value (from tic) captured before method call
%
%   The function never throws.  All errors are captured in the output struct.

    narginchk(4, 4);

    runtime = toc(t_start);   % capture immediately

    % ------------------------------------------------------------------ %
    %  1. Detect hard failure — raw is not a struct                        %
    % ------------------------------------------------------------------ %
    if ~isstruct(raw)
        out = i_failure_output(method_name, cfg, runtime, ...
            'Method did not return a struct.');
        return
    end

    % ------------------------------------------------------------------ %
    %  2. Detect method-reported failure                                   %
    % ------------------------------------------------------------------ %
    if isfield(raw, 'success') && islogical(raw.success) && ~raw.success
        msg = '';
        if isfield(raw, 'message') && ischar(raw.message)
            msg = raw.message;
        end
        out = i_failure_output(method_name, cfg, runtime, ...
            ['Method reported failure. ', msg]);
        out.raw_output = raw;
        return
    end

    % ------------------------------------------------------------------ %
    %  3. Build fields struct from raw                                     %
    % ------------------------------------------------------------------ %
    fields = struct();

    fields.coefficients    = i_safe_get(raw, 'B',            []);
    fields.residuals       = i_safe_get(raw, 'residuals',    []);
    fields.fitted_values   = i_safe_get(raw, 'fitted_values',[]);
    fields.predictions     = i_safe_get(raw, 'predictions',  []);
    fields.selected_terms  = i_build_selected_terms(raw);
    fields.lag_order       = i_safe_get(raw, 'p',            NaN);
    fields.hyperparameters = i_build_hyperparameters(raw, cfg);
    fields.warnings        = i_safe_get_cell(raw, 'warnings');
    fields.diagnostics     = i_safe_get_struct(raw, 'diagnostics');
    fields.message         = i_safe_get(raw, 'message', '');
    fields.raw_output      = raw;

    % Fitted values: derive from data if not supplied by method
    if isempty(fields.fitted_values) && ...
            ~isempty(fields.residuals) && isfield(raw, 'Y0')
        fields.fitted_values = raw.Y0 - fields.residuals;
    end

    % ------------------------------------------------------------------ %
    %  4. Build metadata                                                   %
    % ------------------------------------------------------------------ %
    meta = struct();
    meta.runtime      = runtime;
    meta.T            = i_safe_scalar(raw, 'T');
    meta.K            = i_safe_scalar(raw, 'K');
    meta.T_eff        = i_safe_scalar(raw, 'T_eff');
    meta.toolbox_used = i_safe_get(raw, 'toolbox_used', 'none');
    meta.extra        = i_safe_get_struct(raw, 'meta_extra');

    % Derive T_eff from residuals if not provided
    if isnan(meta.T_eff) && ~isempty(fields.residuals)
        meta.T_eff = size(fields.residuals, 1);
    end

    % ------------------------------------------------------------------ %
    %  5. Assemble canonical output                                        %
    % ------------------------------------------------------------------ %
    out = create_canonical_output(method_name, true, fields, meta);

    % ------------------------------------------------------------------ %
    %  6. Validate and stamp status_code accordingly                       %
    % ------------------------------------------------------------------ %
    [is_valid, report] = validate_canonical_output(out);

    if ~is_valid
        out.success     = false;
        out.status_code = 3;   % 3 = output failed post-hoc validation
        out.message     = strjoin(report.errors, ' | ');
    end

    % Merge any validation warnings
    out.warnings = [out.warnings, report.warnings];
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function out = i_failure_output(method_name, cfg, runtime, msg)
    fields            = struct();
    fields.message    = msg;
    meta              = struct();
    meta.runtime      = runtime;
    meta.T            = NaN;
    meta.K            = NaN;
    meta.T_eff        = NaN;
    meta.toolbox_used = 'none';
    meta.extra        = struct();
    out               = create_canonical_output(method_name, false, fields, meta);
    out.status_code   = 1;
    out.message       = msg;
    if isfield(cfg, 'p'), out.lag_order = cfg.p; end
end

function v = i_safe_get(s, fname, default)
    if isfield(s, fname) && ~isempty(s.(fname))
        v = s.(fname);
    else
        v = default;
    end
end

function v = i_safe_get_cell(s, fname)
    if isfield(s, fname) && iscell(s.(fname))
        v = s.(fname);
    else
        v = {};
    end
end

function v = i_safe_get_struct(s, fname)
    if isfield(s, fname) && isstruct(s.(fname))
        v = s.(fname);
    else
        v = struct();
    end
end

function v = i_safe_scalar(s, fname)
    if isfield(s, fname) && isnumeric(s.(fname)) && isscalar(s.(fname))
        v = s.(fname);
    else
        v = NaN;
    end
end

function terms = i_build_selected_terms(raw)
    % Derive selected_terms from non-zero rows of B if not supplied
    if isfield(raw, 'selected_terms') && iscell(raw.selected_terms)
        terms = raw.selected_terms;
    elseif isfield(raw, 'B') && ~isempty(raw.B)
        nz    = any(raw.B ~= 0, 2);
        idx   = find(nz);
        terms = arrayfun(@(x) sprintf('lag_term_%d', x), idx, ...
                         'UniformOutput', false);
    else
        terms = {};
    end
end

function hp = i_build_hyperparameters(raw, cfg)
    % If wrapper already built a complete hyperparameters struct, use it
    if isfield(raw, 'hyperparameters') && isstruct(raw.hyperparameters) ...
            && ~isempty(fieldnames(raw.hyperparameters))
        hp = raw.hyperparameters;
        return
    end

    % Fallback: scrape individual fields from raw and cfg
    hp = struct();
    for src = {raw, cfg}
        s = src{1};
        for fname = {'lambda', 'alpha', 'rho', 'gamma', 'cv', 'n_folds'}
            f = fname{1};
            if isfield(s, f)
                hp.(f) = s.(f);
            end
        end
    end
end