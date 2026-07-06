% layer1_analysis.m
% Layer 1 Statistical Analysis — Per-configuration paired t-tests.
%
% For each of the 18 experimental configurations, identifies the best
% method (highest mean MCC) and runs paired t-tests between it and every
% other method using the 1000 (or 200 for S5) raw MCC values per method.
% Applies Holm correction within each configuration across the 13 tests.
%
% INPUT:  Merged .mat files in results/task3/ (must contain raw struct)
% OUTPUT: layer1_results/
%           <config_tag>_layer1.csv   — one table per configuration
%           layer1_summary.csv        — all configs stacked in one file
%           layer1_log.txt            — full console log
%
% Author: Dimitris Antonopoulos — AUTH THMMY Thesis
% -------------------------------------------------------------------------

clearvars; clc;

% =========================================================================
%  0. PATHS
% =========================================================================
results_dir = 'C:\Users\dimit\OneDrive\ΗΜΜΥ\Διπλωματική\Chapter 2 - methods software\simulations\results\task3\';
output_dir  = fullfile(results_dir, 'layer1_results');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

log_file = fullfile(output_dir, 'layer1_log.txt');
fid_log  = fopen(log_file, 'w', 'n', 'UTF-8');

function dual_print(fid, fmt, varargin)
    msg = sprintf(fmt, varargin{:});
    fprintf('%s', msg);
    fprintf(fid, '%s', msg);
end

dual_print(fid_log, '=========================================================\n');
dual_print(fid_log, '  LAYER 1 ANALYSIS — Per-Configuration Paired t-tests\n');
dual_print(fid_log, '  Run: %s\n', datestr(now));
dual_print(fid_log, '=========================================================\n\n');

% =========================================================================
%  1. CONFIGURATION LIST
% =========================================================================
configs = {
    'S1_N100_pmax5',  'S1';
    'S1_N100_pmax10', 'S1';
    'S2_N50_pmax5',   'S2';
    'S2_N100_pmax5',  'S2';
    'S2_N1000_pmax5', 'S2';
    'S3_N100_pmax4',  'S3';
    'S3_N100_pmax5',  'S3';
    'S3_N200_pmax4',  'S3';
    'S3_N200_pmax5',  'S3';
    'S3_N500_pmax4',  'S3';
    'S3_N500_pmax5',  'S3';
    'S4_K5_N512',     'S4';
    'S4_K5_N1024',    'S4';
    'S4_K10_N512',    'S4';
    'S4_K10_N1024',   'S4';
    'S4_K20_N512',    'S4';
    'S5_N200_pmax3',  'S5';
    'S5_N500_pmax3',  'S5';
};

n_configs = size(configs, 1);

% Method labels — must match field names in raw struct exactly
method_keys = {
    'Full',        'Full VAR';
    'LASSO',       'LASSO';
    'ElasticNet',  'Elastic Net';
    'mBTS',        'mBTS';
    'AdapLASSO',   'Adaptive LASSO';
    'SCAD',        'SCAD';
    'MCP',         'MCP';
    'HLAG_OO',     'HLAG-OO';
    'HLAG_C',      'HLAG-C';
    'BGR',         'BGR';
    'PDS_LM',      'PDS-LM';
    'msVAR',       'msVAR';
    'gLASSO_V',    'Group LASSO-V';
    'gLASSO_L',    'Group LASSO-L';
};

n_methods   = size(method_keys, 1);
alpha       = 0.05;

% =========================================================================
%  2. SUMMARY TABLE — pre-allocate
% =========================================================================
% Columns: config, dgp, method_label, mean_MCC, delta_vs_best,
%          p_raw, p_holm, sig_worse, is_reference
summary_rows = {};

% =========================================================================
%  3. MAIN LOOP OVER CONFIGURATIONS
% =========================================================================
n_errors = 0;

for ci = 1:n_configs
    tag = configs{ci, 1};
    dgp = configs{ci, 2};

    dual_print(fid_log, '\n---------------------------------------------------------\n');
    dual_print(fid_log, '[%d/%d] Config: %s\n', ci, n_configs, tag);

    mat_file = fullfile(results_dir, [tag '.mat']);
    if ~exist(mat_file, 'file')
        dual_print(fid_log, '  ERROR: File not found — skipping.\n');
        n_errors = n_errors + 1;
        continue
    end

    % Load raw struct
    try
        data = load(mat_file, 'raw');
        raw  = data.raw;
    catch ME
        dual_print(fid_log, '  ERROR loading raw struct: %s\n', ME.message);
        n_errors = n_errors + 1;
        continue
    end

    % ------------------------------------------------------------------
    %  3a. Extract MCC vectors for all methods
    % ------------------------------------------------------------------
    mcc_data   = struct();   % mcc_data.(key) = vector of MCC values
    mean_mcc   = nan(n_methods, 1);
    n_valid    = nan(n_methods, 1);

    for mi = 1:n_methods
        key   = method_keys{mi, 1};
        label = method_keys{mi, 2};

        if ~isfield(raw, key)
            dual_print(fid_log, '  WARNING: Field raw.%s not found.\n', key);
            continue
        end

        vec = raw.(key).MCC;

        % Flatten to column vector
        vec = vec(:);

        % Remove NaN
        vec_clean = vec(~isnan(vec));

        if isempty(vec_clean)
            dual_print(fid_log, '  WARNING: %s has no valid MCC values.\n', label);
            continue
        end

        mcc_data.(key) = vec_clean;
        mean_mcc(mi)   = mean(vec_clean);
        n_valid(mi)    = length(vec_clean);
    end

    % ------------------------------------------------------------------
    %  3b. Identify best method (highest mean MCC, ignoring NaN)
    % ------------------------------------------------------------------
    [best_mean, best_idx] = max(mean_mcc);

    if isnan(best_mean)
        dual_print(fid_log, '  ERROR: No valid methods found — skipping config.\n');
        n_errors = n_errors + 1;
        continue
    end

    best_key   = method_keys{best_idx, 1};
    best_label = method_keys{best_idx, 2};
    best_vec   = mcc_data.(best_key);

    dual_print(fid_log, '  Reference (best): %s  (mean MCC = %.4f,  n = %d)\n', ...
        best_label, best_mean, length(best_vec));

    % ------------------------------------------------------------------
    %  3c. Run paired t-tests: best vs every other method
    %      Only between methods where paired comparison is valid
    %      (same number of realizations — use min overlap length)
    % ------------------------------------------------------------------
    p_raw    = nan(n_methods, 1);   % raw p-values
    t_stat   = nan(n_methods, 1);   % t statistics
    n_paired = nan(n_methods, 1);   % number of paired observations used

    non_ref_indices = [];

    for mi = 1:n_methods
        if mi == best_idx, continue; end

        key = method_keys{mi, 1};
        if ~isfield(mcc_data, key), continue; end

        vec = mcc_data.(key);

        % Use the shorter length for pairing (should always be equal
        % within a config, but guard against any edge case)
        n_pair = min(length(best_vec), length(vec));
        if n_pair < 3
            dual_print(fid_log, '  WARNING: %s has fewer than 3 paired obs — skipping.\n', ...
                method_keys{mi,2});
            continue
        end

        x = best_vec(1:n_pair);
        y = vec(1:n_pair);

        [~, p, ~, stats] = ttest(x, y);   % paired t-test: H1: best > other

        p_raw(mi)    = p;
        t_stat(mi)   = stats.tstat;
        n_paired(mi) = n_pair;
        non_ref_indices(end+1) = mi; %#ok<AGROW>
    end

    % ------------------------------------------------------------------
    %  3d. Holm correction — within this configuration only
    % ------------------------------------------------------------------
    p_holm = nan(n_methods, 1);

    if ~isempty(non_ref_indices)
        raw_p_subset = p_raw(non_ref_indices);
        adj_p        = i_holm(raw_p_subset);

        for k = 1:length(non_ref_indices)
            p_holm(non_ref_indices(k)) = adj_p(k);
        end
    end

    % ------------------------------------------------------------------
    %  3e. Build output table for this configuration
    % ------------------------------------------------------------------
    dual_print(fid_log, '\n  %-18s  %8s  %10s  %10s  %10s  %s\n', ...
        'Method', 'Mean MCC', 'Delta', 'p_raw', 'p_Holm', 'Sig.worse?');
    dual_print(fid_log, '  %s\n', repmat('-', 1, 72));

    % Collect rows sorted by mean MCC descending
    [~, sort_idx] = sort(mean_mcc, 'descend', 'MissingPlacement', 'last');

    csv_rows = {};

    for k = 1:length(sort_idx)
        mi    = sort_idx(k);
        key   = method_keys{mi, 1};
        label = method_keys{mi, 2};

        if isnan(mean_mcc(mi))
            % Method not available in this config
            delta_str  = 'NaN';
            p_raw_str  = 'NaN';
            p_holm_str = 'NaN';
            sig_str    = 'NaN';
            is_ref     = 0;
            sig_worse  = NaN;
        elseif mi == best_idx
            delta_str  = '—';
            p_raw_str  = '—';
            p_holm_str = '—';
            sig_str    = 'No (reference)';
            is_ref     = 1;
            sig_worse  = 0;
        else
            delta      = mean_mcc(mi) - best_mean;
            delta_str  = sprintf('%.4f', delta);

            if isnan(p_raw(mi))
                p_raw_str  = 'NaN';
                p_holm_str = 'NaN';
                sig_str    = 'NaN';
                sig_worse  = NaN;
            else
                p_raw_str = sprintf('%.4f', p_raw(mi));

                if p_holm(mi) < 0.001
                    p_holm_str = '<0.001';
                else
                    p_holm_str = sprintf('%.4f', p_holm(mi));
                end

                sig_worse = double(p_holm(mi) < alpha);
                if sig_worse
                    sig_str = 'Yes';
                else
                    sig_str = 'No';
                end
            end
            is_ref = 0;
        end

        % Print to log
        if mi == best_idx
            dual_print(fid_log, '  %-18s  %8.4f  %10s  %10s  %10s  %s\n', ...
                label, mean_mcc(mi), delta_str, p_raw_str, p_holm_str, sig_str);
        elseif ~isnan(mean_mcc(mi))
            dual_print(fid_log, '  %-18s  %8.4f  %10s  %10s  %10s  %s\n', ...
                label, mean_mcc(mi), delta_str, p_raw_str, p_holm_str, sig_str);
        else
            dual_print(fid_log, '  %-18s  %8s  %10s  %10s  %10s  %s\n', ...
                label, 'NaN', delta_str, p_raw_str, p_holm_str, sig_str);
        end

        % Store for CSV
        if isnan(mean_mcc(mi))
            csv_rows{end+1} = {tag, dgp, label, NaN, NaN, NaN, NaN, NaN, is_ref}; %#ok<AGROW>
        else
            csv_rows{end+1} = {tag, dgp, label, mean_mcc(mi), ...
                mean_mcc(mi) - best_mean, p_raw(mi), p_holm(mi), sig_worse, is_ref}; %#ok<AGROW>
        end

        % Also accumulate in summary
        summary_rows{end+1} = csv_rows{end}; %#ok<AGROW>
    end

    % Write per-config CSV
    csv_path = fullfile(output_dir, [tag '_layer1.csv']);
    fid_csv  = fopen(csv_path, 'w', 'n', 'UTF-8');
    fprintf(fid_csv, 'config,dgp_family,method_label,mean_MCC,delta_vs_best,p_raw,p_holm,sig_worse,is_reference\n');
    for r = 1:length(csv_rows)
        row = csv_rows{r};
        fprintf(fid_csv, '%s,%s,%s,%.6f,%.6f,%s,%s,%s,%d\n', ...
            row{1}, row{2}, row{3}, ...
            i_fmt_num(row{4}), i_fmt_num(row{5}), ...
            i_fmt_p(row{6}), i_fmt_p(row{7}), ...
            i_fmt_sig(row{8}), row{9});
    end
    fclose(fid_csv);
    dual_print(fid_log, '\n  Saved: %s\n', csv_path);

end % end config loop

% =========================================================================
%  4. WRITE SUMMARY CSV
% =========================================================================
% Write summary CSV using writetable for guaranteed UTF-8 encoding
summary_path = fullfile(output_dir, 'layer1_summary.csv');

n_rows = length(summary_rows);
cfg_col    = cell(n_rows, 1);
dgp_col    = cell(n_rows, 1);
meth_col   = cell(n_rows, 1);
mcc_col    = nan(n_rows, 1);
delta_col  = nan(n_rows, 1);
praw_col   = nan(n_rows, 1);
pholm_col  = nan(n_rows, 1);
sigw_col   = nan(n_rows, 1);
isref_col  = zeros(n_rows, 1);

for r = 1:n_rows
    row = summary_rows{r};
    cfg_col{r}   = row{1};
    dgp_col{r}   = row{2};
    meth_col{r}  = row{3};
    if ~isnan(row{4}), mcc_col(r)   = row{4}; end
    if ~isnan(row{5}), delta_col(r) = row{5}; end
    if ~isnan(row{6}), praw_col(r)  = row{6}; end
    if ~isnan(row{7}), pholm_col(r) = row{7}; end
    if ~isnan(row{8}), sigw_col(r)  = row{8}; end
    isref_col(r) = row{9};
end

T = table(cfg_col, dgp_col, meth_col, mcc_col, delta_col, ...
          praw_col, pholm_col, sigw_col, isref_col, ...
    'VariableNames', {'config','dgp_family','method_label','mean_MCC', ...
                      'delta_vs_best','p_raw','p_holm','sig_worse','is_reference'});

writetable(T, summary_path);
dual_print(fid_log, 'Summary CSV saved: %s\n', summary_path);

% =========================================================================
%  5. FINAL SUMMARY REPORT
% =========================================================================
dual_print(fid_log, '\n\n=========================================================\n');
dual_print(fid_log, '  LAYER 1 COMPLETE\n');
dual_print(fid_log, '=========================================================\n');
dual_print(fid_log, '  Configurations processed : %d / %d\n', n_configs - n_errors, n_configs);
dual_print(fid_log, '  Errors / skipped         : %d\n', n_errors);
dual_print(fid_log, '  Per-config CSVs saved to : %s\n', output_dir);
dual_print(fid_log, '  Summary CSV              : %s\n', summary_path);
dual_print(fid_log, '  Log file                 : %s\n', log_file);
dual_print(fid_log, '=========================================================\n');

fclose(fid_log);
fprintf('\nDone. Check layer1_results/ for outputs.\n');


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function p_adj = i_holm(p_raw)
% Holm-Bonferroni correction on a vector of raw p-values.
% Returns adjusted p-values in the same order as input.
    n     = length(p_raw);
    p_adj = nan(size(p_raw));

    % Sort ascending
    [p_sorted, sort_idx] = sort(p_raw);
    p_adj_sorted         = nan(size(p_sorted));

    running_max = 0;
    for k = 1:n
        corrected       = p_sorted(k) * (n - k + 1);
        corrected       = min(corrected, 1);
        running_max     = max(running_max, corrected);
        p_adj_sorted(k) = running_max;
    end

    % Restore original order
    p_adj(sort_idx) = p_adj_sorted;
end


function s = i_fmt_num(x)
% Format a number for CSV — returns 'NaN' string if NaN
    if isnan(x)
        s = 'NaN';
    else
        s = sprintf('%.6f', x);
    end
end


function s = i_fmt_p(x)
% Format p-value for CSV
    if isnan(x)
        s = 'NaN';
    else
        s = sprintf('%.6f', x);
    end
end


function s = i_fmt_sig(x)
% Format significance flag for CSV
    if isnan(x)
        s = 'NaN';
    elseif x == 1
        s = '1';
    else
        s = '0';
    end
end