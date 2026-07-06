function run_mc_batch(config_tag, batch_id, n_batches, n_mc_total)
% RUN_MC_BATCH  Runs a contiguous subset of MC iterations for one config.
%
%   Called by run_sim_task3_hpc.m, which is invoked by the Slurm array job.
%   Each Slurm array task executes one batch (e.g., batch_id=3 out of 20).
%   Partial results are saved as {config_tag}_batch{NNN}.mat and later
%   merged by merge_hpc_results.m into the final .mat file.
%
%   INPUTS
%     config_tag  : string, e.g. 'S1_N100_pmax5'
%     batch_id    : integer in [1, n_batches] — this task's index
%     n_batches   : total number of batches (= Slurm array size)
%     n_mc_total  : total MC runs across all batches (default: 1000)
%
%   OUTPUT
%     Saves: results/task3/partial/{config_tag}_batch{NNN}.mat
%
%   USAGE (from Slurm script)
%     matlab -nodisplay -nosplash -r "run_mc_batch('S1_N100_pmax5',3,20,1000); exit"

    if nargin < 4, n_mc_total = 1000; end

    % ------------------------------------------------------------------ %
    %  1. Determine which MC run indices this batch owns                  %
    % ------------------------------------------------------------------ %
    batch_size = ceil(n_mc_total / n_batches);
    mc_start   = (batch_id - 1) * batch_size + 1;
    mc_end     = min(batch_id * batch_size, n_mc_total);
    mc_indices = mc_start : mc_end;
    n_mc_batch = length(mc_indices);

    fprintf('[HPC] config=%s  batch=%d/%d  mc=%d..%d  (%d runs)\n', ...
        config_tag, batch_id, n_batches, mc_start, mc_end, n_mc_batch);

    % ------------------------------------------------------------------ %
    %  2. Parse config tag → DGP parameters                              %
    % ------------------------------------------------------------------ %
    [dgp_cfg, sim_cfg] = i_parse_config(config_tag, n_mc_total);

    % ------------------------------------------------------------------ %
    %  3. Build output path for this batch                                %
    % ------------------------------------------------------------------ %
    partial_dir = fullfile(sim_cfg.results_dir, 'partial');
    if ~exist(partial_dir, 'dir'), mkdir(partial_dir); end

    batch_tag  = sprintf('%s_batch%03d', config_tag, batch_id);
    save_path  = fullfile(partial_dir, [batch_tag '.mat']);

    if exist(save_path, 'file')
        fprintf('[HPC] Partial file already exists — skipping: %s\n', save_path);
        return
    end

    % ------------------------------------------------------------------ %
    %  4. Add framework paths (Linux)                                     %
    % ------------------------------------------------------------------ %
    i_add_paths(sim_cfg.framework_root);

    % ------------------------------------------------------------------ %
    %  5. Build method configs                                            %
    % ------------------------------------------------------------------ %
    K     = dgp_cfg.K;
    if isfield(dgp_cfg, 'pmax')
        pmax = dgp_cfg.pmax;
    else
        pmax = 5;
    end
    method_cfgs = define_method_configs(pmax, K);
    n_methods   = length(method_cfgs);

    % ------------------------------------------------------------------ %
    %  6. Pre-generate B_true for S3 (fixed system, only noise varies)   %
    % ------------------------------------------------------------------ %
    B_true_s3 = [];
    if strcmp(dgp_cfg.name, 'S3')
        [~, B_true_s3, ~] = generate_S3(dgp_cfg.N, dgp_cfg.seed_system, 1, dgp_cfg.burnin);
    end

    % ------------------------------------------------------------------ %
    %  7. Initialise storage for this batch                               %
    % ------------------------------------------------------------------ %
    metric_names = {'sensitivity','specificity','MCC','FM','HD', ...
                    'TP','FP','TN','FN'};
    raw = struct();
    for m = 1:n_methods
        lbl = method_cfgs(m).label;
        for mn = 1:length(metric_names)
            raw.(lbl).(metric_names{mn}) = NaN(n_mc_batch, 1);
        end
        raw.(lbl).runtime    = NaN(n_mc_batch, 1);
        raw.(lbl).n_failures = 0;
    end

    gc_true = [];

    % ------------------------------------------------------------------ %
    %  8. MC loop (sequential — parallelism comes from Slurm array)      %
    % ------------------------------------------------------------------ %
    t_batch = tic;

    for local_idx = 1:n_mc_batch
        mc = mc_indices(local_idx);   % global MC run index (= random seed)

        % Generate data
        switch dgp_cfg.name
            case 'S1'
                [Y, ~, gc_true] = generate_S1(dgp_cfg.N, mc, dgp_cfg.burnin);
            case 'S2'
                [Y, ~, gc_true] = generate_S2(dgp_cfg.N, mc, dgp_cfg.burnin);
            case 'S3'
                [Y, ~, gc_true] = generate_S3(dgp_cfg.N, dgp_cfg.seed_system, mc, dgp_cfg.burnin);
            case 'S4'
                [Y, gc_true]    = generate_S4(dgp_cfg.N, dgp_cfg.K, dgp_cfg.C, mc, dgp_cfg.burnin);
            case 'S5'
                [Y, ~, gc_true] = generate_S5(dgp_cfg.N, dgp_cfg.seed_system, mc, dgp_cfg.burnin);
        end
        
        % Run all methods
        for m = 1:n_methods
            result = run_single_experiment(Y, gc_true, method_cfgs(m), pmax, sim_cfg.alpha_fdr);
            lbl    = method_cfgs(m).label;

            raw.(lbl).runtime(local_idx) = result.runtime;
            if ~result.success
                raw.(lbl).n_failures = raw.(lbl).n_failures + 1;
                continue
            end
            mf = result.metrics;
            for mn = 1:length(metric_names)
                raw.(lbl).(metric_names{mn})(local_idx) = mf.(metric_names{mn});
            end
        end

        % Progress every 10 runs
        if mod(local_idx, 10) == 0
            elapsed = toc(t_batch);
            eta     = elapsed / local_idx * (n_mc_batch - local_idx);
            fprintf('[HPC] batch %d | run %d/%d | elapsed %.1fs | ETA %.1fs\n', ...
                batch_id, local_idx, n_mc_batch, elapsed, eta);
        end
    end

    % ------------------------------------------------------------------ %
    %  9. Save partial results                                            %
    % ------------------------------------------------------------------ %
    batch_info.config_tag = config_tag;
    batch_info.batch_id   = batch_id;
    batch_info.n_batches  = n_batches;
    batch_info.mc_indices = mc_indices;
    batch_info.n_mc_batch = n_mc_batch;

    save(save_path, 'raw', 'gc_true', 'dgp_cfg', 'sim_cfg', 'batch_info', '-v7.3');
    fprintf('[HPC] Saved: %s  (%.1fs total)\n', save_path, toc(t_batch));
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function [dgp_cfg, sim_cfg] = i_parse_config(tag, n_mc_total)
    % Parse 'S1_N100_pmax5' → dgp_cfg struct
    % Mirrors the logic in run_sim_task3.m

    % Defaults
    sim_cfg.alpha_fdr    = 0.05;
    sim_cfg.burnin       = 200;
    sim_cfg.n_mc_total   = n_mc_total;
    sim_cfg.results_dir  = fullfile(fileparts(mfilename('fullpath')), 'results', 'task3');
    sim_cfg.framework_root = '/home/d/ddantono/svar';

    dgp_cfg.name   = '';
    dgp_cfg.burnin = 200;

    parts = strsplit(tag, '_');

    % DGP name (first token)
    dgp_cfg.name = parts{1};   % S1, S2, S3, S4

    % Parse remaining tokens
    for p = 2:length(parts)
        tok = parts{p};
        if startsWith(tok, 'N')
            dgp_cfg.N = str2double(tok(2:end));
        elseif startsWith(tok, 'pmax')
            dgp_cfg.pmax = str2double(tok(5:end));
        elseif startsWith(tok, 'K')
            dgp_cfg.K = str2double(tok(2:end));
        end
    end

    % DGP-specific defaults
    switch dgp_cfg.name
        case 'S1'
            dgp_cfg.K           = 5;
            dgp_cfg.true_order  = 4;
        case 'S2'
            dgp_cfg.K           = 4;
            dgp_cfg.true_order  = 5;
        case 'S3'
            dgp_cfg.K           = 20;
            dgp_cfg.true_order  = 3;
            dgp_cfg.seed_system = 42;
        case 'S4'
            dgp_cfg.C = 0.5;
            if ~isfield(dgp_cfg, 'K'), dgp_cfg.K = 5; end

        case 'S5'
            dgp_cfg.K           = 50;
            dgp_cfg.true_order  = 2;
            dgp_cfg.seed_system = 99;
    end
end

function i_add_paths(framework_root)
    % Add all required paths (Linux version — no Windows R PATH)
    sims = fullfile(framework_root, 'simulations');
    addpath(fullfile(framework_root, 'core'));
    addpath(fullfile(framework_root, 'utils'));
    addpath(sims);
    addpath(fullfile(sims, 'dgp'));
    addpath(fullfile(sims, 'cgci'));

    % Methods
    methods_dir = fullfile(framework_root, 'methods');
    d = dir(methods_dir);
    for i = 1:length(d)
        if d(i).isdir && ~startsWith(d(i).name, '.')
            addpath(fullfile(methods_dir, d(i).name));
        end
    end
    % R bridge
    addpath(fullfile(framework_root, 'core', 'r_bridge'));
end
