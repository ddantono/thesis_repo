function [is_valid, report] = validate_canonical_output(out)
% VALIDATE_CANONICAL_OUTPUT  Checks a canonical struct for completeness.
%
%   [is_valid, report] = validate_canonical_output(out)
%
%   is_valid : logical scalar
%   report   : struct with fields
%                .passed   — logical scalar
%                .errors   — cell array of error strings
%                .warnings — cell array of warning strings
%                .checks   — struct of per-field pass/fail flags

    report.passed   = true;
    report.errors   = {};
    report.warnings = {};
    report.checks   = struct();

    % ------------------------------------------------------------------ %
    %  1. Top-level type guard                                             %
    % ------------------------------------------------------------------ %
    if ~isstruct(out)
        report.passed = false;
        report.errors{end+1} = 'Output is not a struct.';
        is_valid = false;
        return
    end

    % ------------------------------------------------------------------ %
    %  2. Mandatory field presence                                         %
    % ------------------------------------------------------------------ %
    mandatory = { ...
        'method_name', 'success', 'status_code', 'message', ...
        'coefficients', 'selected_terms', 'lag_order', 'hyperparameters', ...
        'residuals', 'fitted_values', 'predictions', ...
        'runtime', 'warnings', 'diagnostics', 'raw_output', 'metadata'};

    for i = 1:numel(mandatory)
        fname  = mandatory{i};
        present = isfield(out, fname);
        report.checks.(fname) = present;
        if ~present
            report.errors{end+1} = sprintf('Missing mandatory field: "%s".', fname);
            report.passed = false;
        end
    end

    % Abort further checks if fields are missing
    if ~report.passed
        is_valid = false;
        return
    end

    % ------------------------------------------------------------------ %
    %  3. Type checks on critical fields                                   %
    % ------------------------------------------------------------------ %

    % method_name
    if ~ischar(out.method_name) && ~isstring(out.method_name)
        report.errors{end+1} = 'method_name must be char or string.';
        report.passed = false;
    end

    % success
    if ~islogical(out.success) || ~isscalar(out.success)
        report.errors{end+1} = 'success must be a logical scalar.';
        report.passed = false;
    end

    % status_code — must be non-negative integer scalar
    if ~isnumeric(out.status_code) || ~isscalar(out.status_code) ...
            || out.status_code < 0 || out.status_code ~= floor(out.status_code)
        report.errors{end+1} = 'status_code must be a non-negative integer scalar.';
        report.passed = false;
    end

    % success/status_code consistency
    if islogical(out.success) && isnumeric(out.status_code)
        if out.success && out.status_code == 1
            report.errors{end+1} = ...
                'Inconsistency: success=true but status_code=1 (failure).';
            report.passed = false;
        end
        if ~out.success && out.status_code == 0
            report.errors{end+1} = ...
                'Inconsistency: success=false but status_code=0 (ok).';
            report.passed = false;
        end
    end

    % coefficients — numeric or empty
    if ~isnumeric(out.coefficients)
        report.errors{end+1} = 'coefficients must be numeric (or []).';
        report.passed = false;
    end

    % residuals — numeric or empty
    if ~isnumeric(out.residuals)
        report.errors{end+1} = 'residuals must be numeric (or []).';
        report.passed = false;
    end

    % fitted_values — numeric or empty
    if ~isnumeric(out.fitted_values)
        report.errors{end+1} = 'fitted_values must be numeric (or []).';
        report.passed = false;
    end

    % lag_order — scalar positive int or NaN
    ok_p = isnumeric(out.lag_order) && isscalar(out.lag_order) ...
           && (isnan(out.lag_order) || ...
               (out.lag_order >= 1 && out.lag_order == floor(out.lag_order)));
    if ~ok_p
        report.errors{end+1} = ...
            'lag_order must be a positive integer scalar or NaN.';
        report.passed = false;
    end

    % runtime — non-negative scalar or NaN
    ok_rt = isnumeric(out.runtime) && isscalar(out.runtime) ...
            && (isnan(out.runtime) || out.runtime >= 0);
    if ~ok_rt
        report.errors{end+1} = 'runtime must be a non-negative scalar or NaN.';
        report.passed = false;
    end

    % warnings / selected_terms — cell arrays
    if ~iscell(out.warnings)
        report.errors{end+1} = 'warnings must be a cell array.';
        report.passed = false;
    end
    if ~iscell(out.selected_terms)
        report.errors{end+1} = 'selected_terms must be a cell array.';
        report.passed = false;
    end

    % hyperparameters / diagnostics — structs
    if ~isstruct(out.hyperparameters)
        report.errors{end+1} = 'hyperparameters must be a struct.';
        report.passed = false;
    end
    if ~isstruct(out.diagnostics)
        report.errors{end+1} = 'diagnostics must be a struct.';
        report.passed = false;
    end
    
    % Mandatory diagnostics fields
    if isstruct(out.diagnostics) && out.success
        diag_required = {'sparsity', 'n_nonzero', 'mse_per_eq', 'oos_msfe'};
        for i = 1:numel(diag_required)
            if ~isfield(out.diagnostics, diag_required{i})
                report.warnings{end+1} = sprintf( ...
                    'diagnostics.%s is missing.', diag_required{i});
            end
        end
    end
    
    % metadata — struct with sub-fields
    if ~isstruct(out.metadata)
        report.errors{end+1} = 'metadata must be a struct.';
        report.passed = false;
    else
        meta_required = {'method_name', 'created_at', 'T', 'K', 'T_eff'};
        for i = 1:numel(meta_required)
            if ~isfield(out.metadata, meta_required{i})
                report.warnings{end+1} = sprintf( ...
                    'metadata.%s is missing.', meta_required{i});
            end
        end
    end

    % ------------------------------------------------------------------ %
    %  4. Soft checks → warnings only                                      %
    % ------------------------------------------------------------------ %

    % Warn if coefficients present but no selected_terms
    if ~isempty(out.coefficients) && isempty(out.selected_terms)
        report.warnings{end+1} = ...
            'coefficients are present but selected_terms is empty.';
    end

    % Warn if success=true but coefficients empty
    if out.success && isempty(out.coefficients)
        report.warnings{end+1} = ...
            'Method succeeded but coefficients field is empty.';
    end

    % Warn if method_output contains NaN/Inf in residuals
    if ~isempty(out.residuals) && any(~isfinite(out.residuals(:)))
        report.warnings{end+1} = 'residuals contain NaN or Inf values.';
    end

    % ------------------------------------------------------------------ %
    %  5. Final verdict                                                    %
    % ------------------------------------------------------------------ %
    is_valid = report.passed;
end