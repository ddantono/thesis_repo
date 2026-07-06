function patch_msvar_results(config_tag)
% PATCH_MSVAR_RESULTS  Merges partial_msvar files and patches the main .mat.
%
%   Replaces only the msVAR field in the merged .mat file.
%   All other method results are completely untouched.
%
%   Run AFTER all msVAR patch batches complete for this config.
%   The main merged .mat file must already exist (run merge_hpc_results first).
%
%   USAGE:
%     patch_msvar_results('S1_N100_pmax5')

    this_dir    = fileparts(mfilename('fullpath'));   % simulations\msvar_patch\
    sim_dir     = fileparts(this_dir);                % simulations\
    results_dir = fullfile(sim_dir, 'results', 'task3');
    partial_dir = fullfile(results_dir, 'partial_msvar');
    main_file   = fullfile(results_dir, [config_tag '.mat']);

    % 1. Verify main merged file exists
    if ~exist(main_file, 'file')
        error('[msVAR_patch] Main merged file not found: %s\nRun merge_hpc_results first.', ...
            main_file);
    end

    % 2. Find patch partial files
    files = dir(fullfile(partial_dir, sprintf('%s_batch*.mat', config_tag)));
    if isempty(files)
        error('[msVAR_patch] No partial_msvar files found for: %s', config_tag);
    end
    [~, idx] = sort({files.name});
    files = files(idx);
    fprintf('[msVAR_patch] Found %d partial files for %s\n', length(files), config_tag);

    % 3. Concatenate raw metrics
    metric_names = {'sensitivity','specificity','MCC','FM','HD','TP','FP','TN','FN'};

    all_raw.msVAR.runtime    = [];
    all_raw.msVAR.n_failures = 0;
    for mn = 1:length(metric_names)
        all_raw.msVAR.(metric_names{mn}) = [];
    end

    for fi = 1:length(files)
        d = load(fullfile(partial_dir, files(fi).name), 'raw');
        fprintf('[msVAR_patch] Loaded %s\n', files(fi).name);
        for mn = 1:length(metric_names)
            all_raw.msVAR.(metric_names{mn}) = ...
                [all_raw.msVAR.(metric_names{mn}); d.raw.msVAR.(metric_names{mn})];
        end
        all_raw.msVAR.runtime    = [all_raw.msVAR.runtime; d.raw.msVAR.runtime];
        all_raw.msVAR.n_failures = all_raw.msVAR.n_failures + d.raw.msVAR.n_failures;
    end

    % 4. Verify run count
    n_runs = length(all_raw.msVAR.MCC);
    fprintf('[msVAR_patch] Total MC runs: %d\n', n_runs);
    if n_runs ~= 1000 && n_runs ~= 200
        warning('[msVAR_patch] Unexpected run count: %d. Check for missing batches.', ...
            n_runs);
    end

    % 5. Compute aggregated statistics
    new_agg_msvar = struct();
    for mn = 1:length(metric_names)
        vals = all_raw.msVAR.(metric_names{mn});
        new_agg_msvar.(metric_names{mn}).mean = mean(vals, 'omitnan');
        new_agg_msvar.(metric_names{mn}).std  = std(vals,  'omitnan');
        new_agg_msvar.(metric_names{mn}).all  = vals;
    end
    new_agg_msvar.runtime.mean = mean(all_raw.msVAR.runtime, 'omitnan');
    new_agg_msvar.runtime.std  = std(all_raw.msVAR.runtime,  'omitnan');
    new_agg_msvar.n_failures   = all_raw.msVAR.n_failures;

    % 6. Load existing agg, replace only msVAR field, save with -append
    data     = load(main_file, 'agg');
    agg      = data.agg;
    agg.msVAR = new_agg_msvar;
    save(main_file, 'agg', '-append');
    fprintf('[msVAR_patch] Patched and saved: %s\n', main_file);

    % 7. Print updated table
    fprintf('\n--- Patched Results: %s ---\n', config_tag);
    fprintf('%-14s  SENS   SPEC    MCC     FM      HD\n', 'Method');
    fprintf('%s\n', repmat('-',1,58));
    methods = fieldnames(agg);
    for i = 1:length(methods)
        m = methods{i};
        s = agg.(m);
        if isnan(s.MCC.mean)
            fprintf('%-14s  NaN\n', m);
        else
            fprintf('%-14s  %.3f  %.3f  %.3f  %.3f  %.3f\n', ...
                m, s.sensitivity.mean, s.specificity.mean, ...
                s.MCC.mean, s.FM.mean, s.HD.mean);
        end
    end
end