%RUN_SIM_TASK3  Master simulation script for Task 3 (GC estimation study).
%
%  Replicates and extends the simulation study of:
%    Siggiridou & Kugiumtzis (2016). IEEE TSP 64(7), Tables II, III, IV, Fig.4.
%
%  PATHS ARE SET UP AUTOMATICALLY based on this file's location.
%  This script assumes the folder structure:
%
%    sparse_var_framework/           <- framework root (run_framework.m lives here)
%    └── simulations/
%        └── run_sim_task3.m         <- this file
%
%  If your structure is different, set FRAMEWORK_ROOT manually in the
%  USER SETTINGS section below.
%
%  QUICK TEST (no R needed, ~seconds):
%    Set N_MC_OVERRIDE = 10
%    Set CONFIGS_OVERRIDE = {'S1_N100_pmax5'}
%    Run this script — Full baseline will show results, others need R/framework
%
%  FULL VALIDATION (framework + R required, ~minutes per config):
%    Set N_MC_OVERRIDE = 1000
%    Set CONFIGS_OVERRIDE = {'S1_N100_pmax5'}
%    Check mBTS row: expect Sens~0.82, Spec~0.94, MCC~0.78 (Table II, paper)
%
%  FULL RUN (overnight):
%    Set both overrides to []

% =========================================================================
%  USER SETTINGS
% =========================================================================

% Override N_MC for quick testing. Set to [] to use per-config defaults (1000).
N_MC_OVERRIDE = [];

% Override which configs to run. Set to {} to run all.
% Tags: 'S1_N100_pmax5', 'S1_N100_pmax10', 'S2_N50_pmax5', etc.
% For S4: 'S4_K5_N512', 'S4_K10_N1024', etc.
CONFIGS_OVERRIDE = {'S1_N100_pmax10'};

% Significance level for FDR correction (matches paper)
ALPHA_FDR = 0.05;

% Burn-in length for data generation
BURNIN = 500;

% Framework root — set manually if auto-detection fails
FRAMEWORK_ROOT = '';   % leave '' for auto-detection

% =========================================================================
%  AUTOMATIC PATH SETUP
% =========================================================================
sim_dir = fileparts(mfilename('fullpath'));

% Auto-detect framework root (parent of simulations folder)
if isempty(FRAMEWORK_ROOT)
    FRAMEWORK_ROOT = fileparts(sim_dir);
end

% Verify framework root contains run_framework.m
if ~exist(fullfile(FRAMEWORK_ROOT, 'run_framework.m'), 'file')
    error(['run_sim_task3: Cannot find run_framework.m in:\n  %s\n' ...
           'Set FRAMEWORK_ROOT manually in the USER SETTINGS section.'], FRAMEWORK_ROOT);
end

% Add framework paths
addpath(FRAMEWORK_ROOT);
addpath(fullfile(FRAMEWORK_ROOT, 'core'));
addpath(fullfile(FRAMEWORK_ROOT, 'core', 'r_bridge'));
addpath(fullfile(FRAMEWORK_ROOT, 'utils'));

% Add all method subfolders
methods_dir = fullfile(FRAMEWORK_ROOT, 'methods');
if exist(methods_dir, 'dir')
    d = dir(methods_dir);
    for ii = 1:length(d)
        if d(ii).isdir && d(ii).name(1) ~= '.'
            addpath(fullfile(methods_dir, d(ii).name));
        end
    end
end

% Add simulation subfolders
addpath(fullfile(sim_dir, 'dgp'));
addpath(fullfile(sim_dir, 'cgci'));
addpath(sim_dir);

% Windows: add R to PATH so R-based methods can call Rscript
if ispc
    r_path = 'C:\Program Files\R\R-4.5.3\bin';
    if exist(r_path, 'dir')
        setenv('PATH', [r_path ';' getenv('PATH')]);
        fprintf('[setup] Added R to PATH: %s\n', r_path);
    else
        fprintf('[setup] WARNING: R path not found at %s\n', r_path);
        fprintf('[setup] R-based methods will fail. Update r_path in this script.\n');
    end
end

fprintf('[setup] Framework root: %s\n', FRAMEWORK_ROOT);

% Results output folder
RESULTS_DIR = fullfile(sim_dir, 'results', 'task3');
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

% =========================================================================
%  SIMULATION CONFIGURATIONS
%  Each config matches a specific table/figure in Siggiridou 2016.
% =========================================================================
configs = {};

% ---- S1: Table II ----  VAR(4), K=5, N=100
%   pmax=5  (left panel — one larger than true order 4)
%   pmax=10 (right panel — more over-specified)
for pmax_val = [5, 10]
    c = struct('name','S1', 'N',100, 'pmax',pmax_val, 'N_MC',1000);
    configs{end+1} = c;
end

% ---- S2: Table III ----  VAR(5), K=4, pmax=5 (true order)
for N_val = [50, 100, 1000]
    c = struct('name','S2', 'N',N_val, 'pmax',5, 'N_MC',1000);
    configs{end+1} = c;
end

% ---- S3: Table IV ----  sparse VAR(3), K=20
%   pmax=4 and pmax=5 crossed with N=100, 200, 500
for N_val = [100, 200, 500]
    for pmax_val = [4, 5]
        c = struct('name','S3', 'N',N_val, 'pmax',pmax_val, 'N_MC',1000, 'seed_system',42);
        configs{end+1} = c;
    end
end

% ---- S4: Fig. 4 ----  Henon coupled maps, C=0.5
%   K=5,10,20 crossed with N=512, 1024
for K_val = [5, 10, 20]
    for N_val = [512, 1024]
        c = struct('name','S4', 'N',N_val, 'pmax',5, 'K',K_val, 'C',0.5, 'N_MC',1000);
        configs{end+1} = c;
    end
end

% =========================================================================
%  Apply overrides
% =========================================================================
if ~isempty(CONFIGS_OVERRIDE)
    keep = false(1, length(configs));
    for ci = 1:length(configs)
        c   = configs{ci};
        tag = i_config_tag(c);
        if any(strcmp(CONFIGS_OVERRIDE, tag))
            keep(ci) = true;
        end
    end
    configs = configs(keep);
    fprintf('[setup] Running %d selected config(s).\n', sum(keep));
end

if ~isempty(N_MC_OVERRIDE)
    for ci = 1:length(configs)
        configs{ci}.N_MC = N_MC_OVERRIDE;
    end
    fprintf('[setup] N_MC overridden to %d.\n', N_MC_OVERRIDE);
end

% =========================================================================
%  Main loop
% =========================================================================
fprintf('\n=========================================================\n');
fprintf('  TASK 3 SIMULATION STUDY\n');
fprintf('  Total configs to run: %d\n', length(configs));
fprintf('  Results directory:    %s\n', RESULTS_DIR);
fprintf('=========================================================\n\n');

t_global = tic;

for ci = 1:length(configs)
    dgp_cfg = configs{ci};
    tag     = i_config_tag(dgp_cfg);

    % Skip if already saved (resume support)
    save_file = fullfile(RESULTS_DIR, [tag '.mat']);
    if exist(save_file, 'file')
        fprintf('[%d/%d] SKIPPING %s  (result already saved)\n', ci, length(configs), tag);
        continue;
    end

    fprintf('[%d/%d] Starting %s ...\n', ci, length(configs), tag);

    % Determine K
    switch upper(dgp_cfg.name)
        case 'S1', K_cfg = 5;
        case 'S2', K_cfg = 4;
        case 'S3', K_cfg = 20;
        case 'S4', K_cfg = dgp_cfg.K;
        otherwise, error('Unknown DGP: %s', dgp_cfg.name);
    end

    % Build method configurations
    method_cfgs = define_method_configs(dgp_cfg.pmax, K_cfg);

    % Simulation control
    sim_cfg = struct( ...
        'N_MC',      dgp_cfg.N_MC, ...
        'alpha_fdr', ALPHA_FDR, ...
        'burnin',    BURNIN, ...
        'save_path', RESULTS_DIR, ...
        'verbose',   true);

    % Run
    try
        run_mc_simulation(dgp_cfg, method_cfgs, sim_cfg);
        fprintf('[%d/%d] DONE %s\n\n', ci, length(configs), tag);
    catch ME
        fprintf('[%d/%d] ERROR in %s: %s\n\n', ci, length(configs), tag, ME.message);
    end
end

fprintf('\n=========================================================\n');
fprintf('  ALL CONFIGS COMPLETE.  Total time: %.1f min\n', toc(t_global)/60);
fprintf('  Results in: %s\n', RESULTS_DIR);
fprintf('=========================================================\n\n');

% =========================================================================
%  LOCAL HELPER
% =========================================================================
function tag = i_config_tag(c)
if isfield(c, 'K')
    tag = sprintf('%s_K%d_N%d', c.name, c.K, c.N);
else
    tag = sprintf('%s_N%d_pmax%d', c.name, c.N, c.pmax);
end
end