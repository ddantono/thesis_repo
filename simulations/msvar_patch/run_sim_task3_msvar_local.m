% run_sim_task3_msvar_local.m
% Master script — runs msVAR patch for all 16 configs locally on Windows.
%
% USAGE: Open in MATLAB and run. No arguments needed.
% Runs configs sequentially. Each batch of 50 runs saves immediately.
% Safe to interrupt and resume — already-saved batches are skipped.

% =========================================================================
%  PATH SETUP
% =========================================================================
this_dir       = fileparts(mfilename('fullpath'));   % simulations\msvar_patch\
sim_dir        = fileparts(this_dir);                % simulations\
framework_root = fileparts(sim_dir);                 % Chapter 2 - methods software\

% Add framework paths
addpath(framework_root);
addpath(fullfile(framework_root, 'core'));
addpath(fullfile(framework_root, 'core', 'r_bridge'));
addpath(fullfile(framework_root, 'utils'));
addpath(sim_dir);
addpath(fullfile(sim_dir, 'dgp'));
addpath(fullfile(sim_dir, 'cgci'));
addpath(fullfile(sim_dir, 'msvar_patch'));

% Add all method subfolders
methods_dir = fullfile(framework_root, 'methods');
d = dir(methods_dir);
for ii = 1:length(d)
    if d(ii).isdir && d(ii).name(1) ~= '.'
        addpath(fullfile(methods_dir, d(ii).name));
    end
end

% Windows R path
r_path = 'C:\Program Files\R\R-4.5.3\bin';
if exist(r_path, 'dir')
    setenv('PATH', [r_path ';' getenv('PATH')]);
    fprintf('[setup] R added to PATH: %s\n', r_path);
else
    error('[setup] R not found at %s — update r_path in this script.', r_path);
end

fprintf('[setup] Framework root: %s\n', framework_root);

% =========================================================================
%  CONFIG LIST — matches HPC exactly, excluding S4_K20_N1024
% =========================================================================
configs = {
    'S1_N100_pmax5',   ...  % 1
    'S1_N100_pmax10',  ...  % 2
    'S2_N50_pmax5',    ...  % 3
    'S2_N100_pmax5',   ...  % 4
    'S2_N1000_pmax5',  ...  % 5
    'S3_N100_pmax4',   ...  % 6
    'S3_N100_pmax5',   ...  % 7
    'S3_N200_pmax4',   ...  % 8
    'S3_N200_pmax5',   ...  % 9
    'S3_N500_pmax4',   ...  % 10
    'S3_N500_pmax5',   ...  % 11
    'S4_K5_N512',      ...  % 12
    'S4_K5_N1024',     ...  % 13
    'S4_K10_N512',     ...  % 14
    'S4_K10_N1024',    ...  % 15
    'S4_K20_N512',     ...  % 16
    'S5_N200_pmax3',   ...  % 17
    'S5_N500_pmax3'    ...  % 18
};

n_batches  = 20;
n_mc_total = 200;

% =========================================================================
%  MAIN LOOP
% =========================================================================
fprintf('\n=========================================================\n');
fprintf('  msVAR PATCH — LOCAL WINDOWS RUN\n');
fprintf('  Total configs: %d  |  Batches per config: %d\n', ...
    length(configs), n_batches);
fprintf('=========================================================\n\n');

t_global = tic;

for ci = 1:length(configs)
    tag = configs{ci};
    fprintf('\n[%d/%d] Config: %s\n', ci, length(configs), tag);
    t_config = tic;

    for batch_id = 1:n_batches
        try
            run_msvar_patch_batch(tag, batch_id, n_batches, n_mc_total);
        catch ME
            fprintf('[ERROR] config=%s batch=%d: %s\n', tag, batch_id, ME.message);
        end
    end

    fprintf('[%d/%d] Config %s complete (%.1f min)\n', ...
        ci, length(configs), tag, toc(t_config)/60);
end

fprintf('\n=========================================================\n');
fprintf('  ALL CONFIGS COMPLETE. Total time: %.1f min\n', toc(t_global)/60);
fprintf('=========================================================\n\n');
fprintf('Next step: run patch_msvar_results for each config.\n');