function out = create_canonical_output(method_name, success, fields, meta)
% CREATE_CANONICAL_OUTPUT  Assembles the canonical result struct.
%
%   out = create_canonical_output(method_name, success, fields, meta)
%
%   method_name : char/string — method identifier
%   success     : logical scalar
%   fields      : struct — populated fields from the method wrapper
%   meta        : struct — runtime, data dimensions, hyperparameters, etc.
%
%   Missing or empty fields are filled with canonical sentinel values.
%   All dimension checks are performed here before the struct is sealed.

    % ------------------------------------------------------------------ %
    %  0. Argument guards                                                  %
    % ------------------------------------------------------------------ %
    narginchk(4, 4);

    assert(ischar(method_name) || isstring(method_name), ...
        '[create_canonical_output] method_name must be char or string.');
    method_name = char(method_name);

    assert(islogical(success) && isscalar(success), ...
        '[create_canonical_output] success must be a logical scalar.');

    assert(isstruct(fields), ...
        '[create_canonical_output] fields must be a struct.');

    assert(isstruct(meta), ...
        '[create_canonical_output] meta must be a struct.');

    % ------------------------------------------------------------------ %
    %  1. Default-initialise every mandatory field                         %
    %     Rule: if a field is missing or empty → sentinel value            %
    %       numeric arrays  → []                                           %
    %       logical         → false                                        %
    %       char            → ''                                           %
    %       cell            → {}                                           %
    %       struct          → struct()                                     %
    % ------------------------------------------------------------------ %
    out.method_name    = method_name;
    out.success        = success;
    out.status_code    = i_status_code(success);   % 0 = ok, see helper
    out.message        = i_pull_char(fields,  'message',   '');

    % --- model estimates ------------------------------------------------
    out.coefficients   = i_pull_num(fields,  'coefficients');
    out.selected_terms = i_pull_cell(fields, 'selected_terms');
    out.lag_order      = i_pull_scalar_int(fields, 'lag_order');
    out.hyperparameters= i_pull_struct(fields,'hyperparameters');

    % --- fit & prediction -----------------------------------------------
    out.residuals      = i_pull_num(fields,  'residuals');
    out.fitted_values  = i_pull_num(fields,  'fitted_values');
    out.predictions    = i_pull_num(fields,  'predictions');

    % --- diagnostics & warnings -----------------------------------------
    out.warnings       = i_pull_cell(fields, 'warnings');
    out.diagnostics    = i_pull_struct(fields,'diagnostics');

    % --- pass-through of unprocessed method output ----------------------
    out.raw_output     = i_pull_any(fields,  'raw_output');

    % --- runtime (seconds) ----------------------------------------------
    out.runtime        = i_pull_nonneg_scalar(meta, 'runtime');

    % --- metadata block -------------------------------------------------
    out.metadata       = i_build_metadata(meta, method_name);

    % ------------------------------------------------------------------ %
    %  2. Cross-field dimension consistency checks                         %
    %     Rule: warn but do NOT error — just flag inconsistency           %
    % ------------------------------------------------------------------ %
    out = i_check_dimensions(out);

end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function v = i_pull_num(s, fname)
    % Returns numeric array or [] sentinel
    if isfield(s, fname) && isnumeric(s.(fname)) && ~isempty(s.(fname))
        v = double(s.(fname));
    else
        v = [];
    end
end

function v = i_pull_scalar_int(s, fname)
    % Returns positive integer scalar or NaN sentinel
    if isfield(s, fname) && isnumeric(s.(fname)) ...
            && isscalar(s.(fname)) && s.(fname) >= 0
        v = round(double(s.(fname)));
    else
        v = NaN;
    end
end

function v = i_pull_nonneg_scalar(s, fname)
    % Returns non-negative scalar double or NaN
    if isfield(s, fname) && isnumeric(s.(fname)) ...
            && isscalar(s.(fname)) && isfinite(s.(fname)) ...
            && s.(fname) >= 0
        v = double(s.(fname));
    else
        v = NaN;
    end
end

function v = i_pull_char(s, fname, default)
    if isfield(s, fname) ...
            && (ischar(s.(fname)) || isstring(s.(fname))) ...
            && ~isempty(s.(fname))
        v = char(s.(fname));
    else
        v = default;
    end
end

function v = i_pull_cell(s, fname)
    if isfield(s, fname) && iscell(s.(fname))
        v = s.(fname);
    else
        v = {};
    end
end

function v = i_pull_struct(s, fname)
    if isfield(s, fname) && isstruct(s.(fname))
        v = s.(fname);
    else
        v = struct();
    end
end

function v = i_pull_any(s, fname)
    % Pass-through — no type constraint for raw_output
    if isfield(s, fname)
        v = s.(fname);
    else
        v = [];
    end
end

function code = i_status_code(success)
    % 0  — success
    % 1  — failure (method returned success=false)
    % Codes >=2 are set externally (e.g. validate_canonical_output)
    if success
        code = 0;
    else
        code = 1;
    end
end

function meta_out = i_build_metadata(meta, method_name)
    meta_out.method_name  = method_name;
    meta_out.created_at   = datestr(now, 'yyyy-mm-dd HH:MM:SS');  %#ok<TNOW1,DATST>
    meta_out.T            = i_pull_scalar_int(meta, 'T');
    meta_out.K            = i_pull_scalar_int(meta, 'K');
    meta_out.T_eff        = i_pull_scalar_int(meta, 'T_eff');
    meta_out.toolbox_used = i_pull_char(meta, 'toolbox_used', 'none');
    meta_out.extra        = i_pull_struct(meta, 'extra');
end

function out = i_check_dimensions(out)
    % Accumulate dimension warnings; do not throw.
    dim_warn = {};

    T_eff = out.metadata.T_eff;
    K     = out.metadata.K;

    % residuals must be [T_eff x K] if both are known
    if ~isempty(out.residuals) && ~isnan(T_eff) && ~isnan(K)
        [r, c] = size(out.residuals);
        if r ~= T_eff || c ~= K
            dim_warn{end+1} = sprintf( ...
                'residuals size [%d x %d] does not match expected [%d x %d].', ...
                r, c, T_eff, K);
        end
    end

    % fitted_values must match residuals if both present
    if ~isempty(out.fitted_values) && ~isempty(out.residuals)
        if ~isequal(size(out.fitted_values), size(out.residuals))
            dim_warn{end+1} = ...
                'fitted_values and residuals have incompatible sizes.';
        end
    end

    % coefficients: expect [K*p x K]; only check if lag_order is known
    p = out.lag_order;
    if ~isempty(out.coefficients) && ~isnan(p) && ~isnan(K)
        [rb, cb] = size(out.coefficients);
        if rb ~= K*p || cb ~= K
            dim_warn{end+1} = sprintf( ...
                'coefficients size [%d x %d] does not match expected [%d x %d].', ...
                rb, cb, K*p, K);
        end
    end

    % Merge new warnings into any already present
    out.warnings = [out.warnings, dim_warn];

    % If dimension issues found, escalate status_code
    if ~isempty(dim_warn) && out.status_code == 0
        out.status_code = 2;   % 2 = success with warnings
    end
end