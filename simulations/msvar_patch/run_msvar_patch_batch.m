function run_msvar_patch_batch(config_tag, batch_id, n_batches, n_mc_total)
% RUN_MSVAR_PATCH_BATCH  Runs msVAR only for one config/batch locally.
%
%   Targeted patch for msVAR which could not run on HPC due to missing
%   gmp-devel system library. Uses identical seeds to HPC jobs so results
%   are directly substitutable into merged .mat files.
%
%   Saves to simulations/results/task3/partial_msvar/
%
%   USAGE (from run_sim_task3_msvar_local.m):
%     run_msvar_patch_batch('S1_N100_pmax5', 1, 20, 1000)

    if nargin < 4, n_mc_total = 1000; end

    % ------------------------------------------------------------------ %
    %  1. Determine MC indices for this batch                             %
    % ------------------------------------------------------------------ %
    batch_size = ceil(n_mc_total / n_batches);
    mc_start   = (batch_id - 1) * batch_size + 1;
    mc_end     = min(batch_id * batch_size, n_mc_total);
    mc_indices = mc_start : mc_end;
    n_mc_batch = length(mc_indices);

    fprintf('[msVAR_PATCH] config=%s  batch=%d/%d  mc=%d..%d  (%d runs)\n', ...
        config_tag, batch_id, n_batches, mc_start, mc_end, n_mc_batch);

    % ------------------------------------------------------------------ %
    %  2. Parse config tag                                                %
    % ------------------------------------------------------------------ %
    [dgp_cfg, sim_cfg] = i_parse_config(config_tag, n_mc_total);

    % ------------------------------------------------------------------ %
    %  3. Output path                                                     %
    % ------------------------------------------------------------------ %
    partial_dir = fullfile(sim_cfg.results_dir, 'partial_msvar');
    if ~exist(partial_dir, 'dir'), mkdir(partial_dir); end

    save_path = fullfile(partial_dir, ...
        sprintf('%s_batch%03d.mat', config_tag, batch_id));

    if exist(save_path, 'file')
        fprintf('[msVAR_PATCH] Already exists — skipping: %s\n', save_path);
        return
    end

    % ------------------------------------------------------------------ %
    %  4. msVAR config — identical to define_method_configs entry 12     %
    % ------------------------------------------------------------------ %
    if isfield(dgp_cfg, 'pmax')
        pmax = dgp_cfg.pmax;
    else
        pmax = 5;
    end

    msvar_cfg.method = 'msvar';
    msvar_cfg.label  = 'msVAR';
    msvar_cfg.cfg    = struct('method', 'msvar', 'p_seq', 1:pmax);

    % ------------------------------------------------------------------ %
    %  5. Initialise storage                                              %
    % ------------------------------------------------------------------ %
    metric_names = {'sensitivity','specificity','MCC','FM','HD', ...
                    'TP','FP','TN','FN'};

    raw.msVAR.runtime    = NaN(n_mc_batch, 1);
    raw.msVAR.n_failures = 0;
    for mn = 1:length(metric_names)
        raw.msVAR.(metric_names{mn}) = NaN(n_mc_batch, 1);
    end

    % ------------------------------------------------------------------ %
    %  6. Pre-generate S3 topology if needed                             %
    % ------------------------------------------------------------------ %
    if strcmp(dgp_cfg.name, 'S3')
        [~, ~, ~] = generate_S3(dgp_cfg.N, dgp_cfg.seed_system, 1, dgp_cfg.burnin);
    end

    gc_true = [];
    t_batch = tic;

    % ------------------------------------------------------------------ %
    %  7. MC loop — mirrors run_mc_batch.m exactly                       %
    % ------------------------------------------------------------------ %
    for local_idx = 1:n_mc_batch
        mc = mc_indices(local_idx);

        % Generate data using IDENTICAL approach as HPC run_mc_batch.m
        switch dgp_cfg.name
            case 'S1'
                [Y, ~, gc_true] = generate_S1(dgp_cfg.N, mc, dgp_cfg.burnin);
            case 'S2'
                [Y, ~, gc_true] = generate_S2(dgp_cfg.N, mc, dgp_cfg.burnin);
            case 'S3'
                [Y, ~, gc_true] = generate_S3(dgp_cfg.N, dgp_cfg.seed_system, ...
                                               mc, dgp_cfg.burnin);
            case 'S4'
                [Y, gc_true]    = generate_S4(dgp_cfg.N, dgp_cfg.K, ...
                                               dgp_cfg.C, mc, dgp_cfg.burnin);
            case 'S5'
                [Y, ~, gc_true] = generate_S5(dgp_cfg.N, dgp_cfg.seed_system, ...
                                               mc, dgp_cfg.burnin);
        end

        % Run msVAR only
        result = run_single_experiment(Y, gc_true, msvar_cfg, pmax, sim_cfg.alpha_fdr);
        raw.msVAR.runtime(local_idx) = result.runtime;

        if ~result.success
            raw.msVAR.n_failures = raw.msVAR.n_failures + 1;
            fprintf('[msVAR_PATCH] mc=%d FAILED: %s\n', mc, result.msg);
            continue
        end

        for mn = 1:length(metric_names)
            raw.msVAR.(metric_names{mn})(local_idx) = ...
                result.metrics.(metric_names{mn});
        end

        % Progress every 10 runs
        if mod(local_idx, 10) == 0
            elapsed = toc(t_batch);
            eta     = elapsed / local_idx * (n_mc_batch - local_idx);
            fprintf('[msVAR_PATCH] batch %d | run %d/%d | elapsed %.1fs | ETA %.1fs\n', ...
                batch_id, local_idx, n_mc_batch, elapsed, eta);
        end
    end

    % ------------------------------------------------------------------ %
    %  8. Save                                                            %
    % ------------------------------------------------------------------ %
    batch_info.config_tag = config_tag;
    batch_info.batch_id   = batch_id;
    batch_info.n_batches  = n_batches;
    batch_info.mc_indices = mc_indices;
    batch_info.n_mc_batch = n_mc_batch;

    save(save_path, 'raw', 'gc_true', 'dgp_cfg', 'sim_cfg', 'batch_info', '-v7.3');
    fprintf('[msVAR_PATCH] Saved: %s  (%.1fs total)\n', save_path, toc(t_batch));
end


% ======================================================================= %
%  Private helpers                                                         %
% ======================================================================= %

function [dgp_cfg, sim_cfg] = i_parse_config(tag, n_mc_total)
    % Mirrors i_parse_config in run_mc_batch.m exactly

    % Locate simulations folder relative to this file
    this_dir = fileparts(mfilename('fullpath'));  % simulations\msvar_patch\
    sim_dir  = fileparts(this_dir);               % simulations\

    sim_cfg.alpha_fdr    = 0.05;
    sim_cfg.burnin       = 200;
    sim_cfg.n_mc_total   = n_mc_total;
    sim_cfg.results_dir  = fullfile(sim_dir, 'results', 'task3');

    dgp_cfg.name   = '';
    dgp_cfg.burnin = 200;

    parts = strsplit(tag, '_');
    dgp_cfg.name = parts{1};

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

    switch dgp_cfg.name
        case 'S1'
            dgp_cfg.K          = 5;
            dgp_cfg.true_order = 4;
        case 'S2'
            dgp_cfg.K          = 4;
            dgp_cfg.true_order = 5;
        case 'S3'
            dgp_cfg.K           = 20;
            dgp_cfg.true_order  = 3;
            dgp_cfg.seed_system = 42;
        case 'S4'
            dgp_cfg.C = 0.5;
            if ~isfield(dgp_cfg, 'K'), dgp_cfg.K = 5; 
            end
        case 'S5'
            dgp_cfg.K           = 50;
            dgp_cfg.true_order  = 2;
            dgp_cfg.seed_system = 99;   % default seed_system for S5
    end
end