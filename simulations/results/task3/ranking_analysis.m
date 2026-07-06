%% ranking_analysis.m
% Extracts per-realization MCC vectors from merged .mat files and exports
% a clean data structure for the Friedman + Wilcoxon ranking analysis in R.
%
% OUTPUT FILES (saved in results_dir):
%   mcc_matrix.mat  — MATLAB workspace with all data
%   mcc_long.csv    — Long-format CSV for R (config, method, realization, MCC)
%
% Author: Dimitris Dantonopoulos — AUTH THMMY Thesis
% -------------------------------------------------------------------------

clear; clc;

%% ---- 1. PATHS & CONFIGURATION ------------------------------------------

results_dir = 'C:\Users\dimit\OneDrive\ΗΜΜΥ\Διπλωματική\Chapter 2 - methods software\simulations\results\task3\';

% The 16 main configurations (1000 realizations each).
% S5 is treated separately (200 realizations).
configs_main = {
    'S1_N100_pmax5',   'S1_N100_pmax10', ...
    'S2_N50_pmax5',    'S2_N100_pmax5',  'S2_N1000_pmax5', ...
    'S3_N100_pmax4',   'S3_N100_pmax5',  'S3_N200_pmax4', ...
    'S3_N200_pmax5',   'S3_N500_pmax4',  'S3_N500_pmax5', ...
    'S4_K5_N512',      'S4_K5_N1024',    'S4_K10_N512', ...
    'S4_K10_N1024',    'S4_K20_N512'
};

configs_s5 = {'S5_N200_pmax3', 'S5_N500_pmax3'};

% DGP family labels for grouping in plots
dgp_family = {
    'S1','S1', ...
    'S2','S2','S2', ...
    'S3','S3','S3','S3','S3','S3', ...
    'S4','S4','S4','S4','S4'
};

% Method names as stored in .mat files (order must match fieldnames)
method_keys = {'Full','LASSO','ElasticNet','mBTS','AdapLASSO','SCAD', ...
               'MCP','HLAG_OO','HLAG_C','BGR','PDS_LM','msVAR', ...
               'gLASSO_V','gLASSO_L'};

% Display names for plots and tables
method_labels = {'Full VAR','LASSO','Elastic Net','mBTS','Adaptive LASSO', ...
                 'SCAD','MCP','HLAG-OO','HLAG-C','BGR','PDS-LM','msVAR', ...
                 'Group LASSO-V','Group LASSO-L'};

N_main    = numel(configs_main);   % 16
N_s5      = numel(configs_s5);     % 2
N_methods = numel(method_keys);    % 14
N_real    = 1000;
N_real_s5 = 200;

%% ---- 2. LOAD MAIN CONFIGURATIONS (1000 realizations) -------------------

fprintf('Loading main configurations...\n');

% MCC_main: [N_real x N_methods x N_configs]
MCC_main = NaN(N_real, N_methods, N_main);

for c = 1:N_main
    fname = fullfile(results_dir, [configs_main{c} '.mat']);
    fprintf('  Loading %s ... ', configs_main{c});
    d = load(fname, 'raw');
    for m = 1:N_methods
        key = method_keys{m};
        if isfield(d.raw, key) && isfield(d.raw.(key), 'MCC')
            vec = d.raw.(key).MCC(:);
            MCC_main(1:numel(vec), m, c) = vec;
        else
            warning('Missing MCC for method %s in config %s', key, configs_main{c});
        end
    end
    fprintf('OK\n');
end

%% ---- 3. LOAD S5 CONFIGURATIONS (200 realizations) ----------------------

fprintf('Loading S5 configurations...\n');

MCC_s5 = NaN(N_real_s5, N_methods, N_s5);

for c = 1:N_s5
    fname = fullfile(results_dir, [configs_s5{c} '.mat']);
    fprintf('  Loading %s ... ', configs_s5{c});
    d = load(fname, 'raw');
    for m = 1:N_methods
        key = method_keys{m};
        if isfield(d.raw, key) && isfield(d.raw.(key), 'MCC')
            vec = d.raw.(key).MCC(:);
            MCC_s5(1:numel(vec), m, c) = vec;
        else
            warning('Missing MCC for method %s in config %s', key, configs_s5{c});
        end
    end
    fprintf('OK\n');
end

%% ---- 4. COMPUTE MEAN MCC TABLE (for quick inspection) ------------------

fprintf('\n--- Mean MCC per method per config (main) ---\n');
fprintf('%-16s', 'Method');
for c = 1:N_main
    fprintf('  %-18s', configs_main{c});
end
fprintf('\n');

MCC_mean_main = squeeze(nanmean(MCC_main, 1));  % [N_methods x N_configs]

for m = 1:N_methods
    fprintf('%-16s', method_labels{m});
    for c = 1:N_main
        fprintf('  %-18.4f', MCC_mean_main(m, c));
    end
    fprintf('\n');
end

%% ---- 5. EXPORT LONG-FORMAT CSV FOR R -----------------------------------

fprintf('\nExporting long-format CSV for R...\n');

csv_path = fullfile(results_dir, 'mcc_long.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, 'config,dgp_family,method_key,method_label,realization,MCC\n');

% Main configs
for c = 1:N_main
    for m = 1:N_methods
        for r = 1:N_real
            val = MCC_main(r, m, c);
            if ~isnan(val)
                fprintf(fid, '%s,%s,%s,"%s",%d,%.6f\n', ...
                    configs_main{c}, dgp_family{c}, ...
                    method_keys{m}, method_labels{m}, r, val);
            end
        end
    end
end

% S5 configs
for c = 1:N_s5
    for m = 1:N_methods
        for r = 1:N_real_s5
            val = MCC_s5(r, m, c);
            if ~isnan(val)
                fprintf(fid, '%s,S5,%s,"%s",%d,%.6f\n', ...
                    configs_s5{c}, ...
                    method_keys{m}, method_labels{m}, r, val);
            end
        end
    end
end

fclose(fid);
fprintf('  Saved: %s\n', csv_path);

%% ---- 6. SAVE MATLAB WORKSPACE ------------------------------------------

mat_path = fullfile(results_dir, 'mcc_matrix.mat');
save(mat_path, 'MCC_main', 'MCC_s5', 'MCC_mean_main', ...
     'configs_main', 'configs_s5', 'dgp_family', ...
     'method_keys', 'method_labels', 'N_main', 'N_s5', ...
     'N_methods', 'N_real', 'N_real_s5');
fprintf('  Saved: %s\n', mat_path);

fprintf('\nDone. Ready for R ranking analysis.\n');