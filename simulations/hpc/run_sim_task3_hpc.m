% run_sim_task3_hpc.m
%
% HPC MASTER SCRIPT — AUTH Aristotelis Cluster
% =============================================
% This script is called by the Slurm job array submission script.
% It is NOT intended for interactive use on Windows.
%
% HOW IT WORKS:
%   - Slurm launches one MATLAB process per array task
%   - Each task receives its batch_id via the HPC_BATCH_ID environment variable
%   - All tasks run independently and save partial results
%   - After all tasks complete, run merge_hpc_results.m to combine them
%
% WORKFLOW:
%   1. Edit USER SETTINGS below (config list, paths, n_batches)
%   2. Transfer files to HPC (scp or git)
%   3. Run: sbatch submit_task3.sh
%   4. Monitor: squeue -u <username>
%   5. After completion: matlab -nodisplay -r "merge_hpc_results; exit"
%
% IMPORTANT NOTES:
%   - R must be loaded before MATLAB: `module load R/4.x.x` in the Slurm script
%   - Set FRAMEWORK_ROOT below to the absolute path on the HPC filesystem
%   - Rscript must be in PATH (guaranteed if R module is loaded)

% ======================================================================= %
%  USER SETTINGS — edit these before submitting                           %
% ======================================================================= %

% Absolute path to framework root on the HPC filesystem
FRAMEWORK_ROOT = '/home/d/ddantono/svar';

% Total MC runs (1000 for final results, 20 for quick test)
N_MC_TOTAL = 200;

% Number of Slurm array tasks (must match #SBATCH --array in .sh script)
N_BATCHES = 20;   % 20 tasks × 50 runs = 1000 total

% Configs to run — comment out completed ones to skip them
CONFIGS = {
    'S1_N100_pmax5',
    'S1_N100_pmax10',
    'S2_N50_pmax5',
    'S2_N100_pmax5',
    'S2_N1000_pmax5',
    'S3_N100_pmax4',
    'S3_N100_pmax5',
    'S3_N200_pmax4',
    'S3_N200_pmax5',
    'S3_N500_pmax4',
    'S3_N500_pmax5',
    'S4_K5_N512',
    'S4_K5_N1024',
    'S4_K10_N512',
    'S4_K10_N1024',
    'S4_K20_N512',
    'S5_N200_pmax3',
    'S5_N500_pmax3',
};

% ======================================================================= %
%  RUNTIME: get batch_id from environment variable set by Slurm           %
% ======================================================================= %

batch_id_str = getenv('HPC_BATCH_ID');
if isempty(batch_id_str)
    % Fallback: check if passed as MATLAB variable (for testing)
    if ~exist('HPC_BATCH_ID', 'var')
        error(['HPC_BATCH_ID not set. ' ...
               'On HPC: set via Slurm script. ' ...
               'For testing: set HPC_BATCH_ID=1 before running.']);
    end
    batch_id = HPC_BATCH_ID;
else
    batch_id = str2double(batch_id_str);
end

% Config index: which config does this batch belong to?
% Strategy: each config gets its own submission. Set CONFIG_IDX in Slurm script.
config_idx_str = getenv('HPC_CONFIG_IDX');
if isempty(config_idx_str)
    if ~exist('HPC_CONFIG_IDX', 'var')
        config_idx = 1;   % default to first config for testing
        warning('HPC_CONFIG_IDX not set — defaulting to config 1.');
    else
        config_idx = HPC_CONFIG_IDX;
    end
else
    config_idx = str2double(config_idx_str);
end

config_tag = CONFIGS{config_idx};

fprintf('=====================================================\n');
fprintf(' AUTH Aristotelis HPC — Task 3 Simulation\n');
fprintf(' Config [%d/%d]: %s\n', config_idx, length(CONFIGS), config_tag);
fprintf(' Batch:  %d / %d\n', batch_id, N_BATCHES);
fprintf(' MC runs per batch: %d\n', ceil(N_MC_TOTAL / N_BATCHES));
fprintf('=====================================================\n');

% ======================================================================= %
%  Add simulation layer to path                                            %
% ======================================================================= %
addpath(fullfile(FRAMEWORK_ROOT, 'simulations'));
addpath(fullfile(FRAMEWORK_ROOT, 'core'));
% Additional paths added inside run_mc_batch via i_add_paths()

% ======================================================================= %
%  Run this batch                                                          %
% ======================================================================= %
setenv('R_LIBS_USER', '/home/d/ddantono/R/x86_64-pc-linux-gnu-library/4.4');
setenv('PATH', ['/mnt/apps/aristotle/site/linux-rocky9-x86_64/gcc-14.2.0/r-4.4.1-542fb4ll3xdxfahqxewr7eyqnrpwzqsg/bin:', getenv('PATH')]);
run_mc_batch(config_tag, batch_id, N_BATCHES, N_MC_TOTAL);

fprintf('[HPC] Batch %d of config %s complete.\n', batch_id, config_tag);
