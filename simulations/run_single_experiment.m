function result = run_single_experiment(Y, gc_true, method_cfg, pmax, alpha_fdr)
%RUN_SINGLE_EXPERIMENT  Run one sVAR method on one dataset and evaluate GC.
%
%  result = run_single_experiment(Y, gc_true, method_cfg, pmax, alpha_fdr)
%
%  PIPELINE:
%    1. Run sVAR method  ->  canonical output (B_hat, support)
%    2. Extract support  ->  binary [K*pmax x K] mask
%    3. OLS-CGCI         ->  CGCI values + p-values
%    4. FDR correction   ->  gc_est binary matrix
%    5. GC metrics       ->  Sensitivity, Specificity, MCC, FM, HD
%
%  CALLING CONVENTION:
%    This function calls dispatch_method() rather than run_framework().
%    run_framework() is a void orchestrator (saves to disk but returns nothing).
%    dispatch_method() returns the canonical output struct directly.
%
%  INPUTS
%    Y           [T x K]       Multivariate time series.
%    gc_true     [K x K]       True binary GC matrix.
%    method_cfg  struct        Fields: .method, .label, .cfg
%    pmax        positive int  VAR lag order for CGCI computation.
%    alpha_fdr   scalar        FDR significance level (default: 0.05).
%
%  OUTPUT
%    result  struct:
%      .method     char      method label
%      .success    logical
%      .runtime    scalar    method runtime in seconds (NaN if failed)
%      .support    logical   [K*p x K] support matrix
%      .CGCI       [K x K]   CGCI values
%      .pval       [K x K]   p-values (diagonal = NaN)
%      .gc_est     [K x K]   estimated binary GC matrix after FDR
%      .metrics    struct    Sensitivity, Specificity, MCC, FM, HD, TP/FP/TN/FN
%      .msg        char      error message (empty if success)

if nargin < 5 || isempty(alpha_fdr)
    alpha_fdr = 0.05;
end

[~, K] = size(Y);

% Initialise with safe defaults
result.method  = method_cfg.label;
result.success = false;
result.runtime = NaN;
result.support = false(K*pmax, K);
result.CGCI    = zeros(K);
result.pval    = nan(K);
result.gc_est  = zeros(K);
result.metrics = [];
result.msg     = '';

t0 = tic;

try
    % =================================================================
    %  Step 1: Run method -> get canonical output
    % =================================================================

    if strcmp(method_cfg.method, 'full_var')
        % Full VAR baseline: all K*pmax regressors in every equation
        out = struct();
        out.coefficients = ones(K*pmax, K);   % dense = all-ones support
        out.success      = true;
        out.message      = '';
        out.runtime      = 0;

    else
        % Use dispatch_method which returns the canonical output directly.
        % (run_framework is a void orchestrator that saves to disk only.)
        out = dispatch_method(Y, method_cfg.cfg);
        result.runtime = toc(t0);

        if ~out.success
            result.msg = sprintf('dispatch_method failed: %s', out.message);
            return;
        end

        % Ensure B has K*pmax rows (pad with zeros if method used fewer lags)
        if size(out.coefficients, 1) ~= K * pmax
            B_adj = zeros(K*pmax, K);
            n_use = min(size(out.coefficients,1), K*pmax);
            B_adj(1:n_use, :) = out.coefficients(1:n_use, :);
            out.coefficients = B_adj;
        end
    end

    % =================================================================
    %  Step 2: Extract support
    % =================================================================
    support = extract_support(out);

    % =================================================================
    %  Step 3: OLS-CGCI
    % =================================================================
    [CGCI, pval] = compute_cgci_ols(Y, support, pmax);

    % =================================================================
    %  Step 4: FDR correction
    % =================================================================
    gc_est = apply_fdr(pval, alpha_fdr);

    % =================================================================
    %  Step 5: GC metrics
    % =================================================================
    metrics = compute_gc_metrics(gc_est, gc_true);

    % Store
    result.success = true;
    result.runtime = toc(t0);
    result.support = support;
    result.CGCI    = CGCI;
    result.pval    = pval;
    result.gc_est  = gc_est;
    result.metrics = metrics;

catch ME
    result.msg     = ME.message;
    result.runtime = toc(t0);
end
end