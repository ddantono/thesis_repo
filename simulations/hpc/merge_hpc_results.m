function merge_hpc_results(config_tag, delete_partials)
% MERGE_HPC_RESULTS  Combine partial batch .mat files into the final result.
%
%   merge_hpc_results(config_tag)
%   merge_hpc_results(config_tag, true)   % also delete partial files
%
%   Run this AFTER all Slurm array tasks for a config have completed.
%   Produces a .mat file identical in format to run_mc_simulation.m output.
%
%   USAGE
%     % On HPC or locally after downloading partial files:
%     merge_hpc_results('S1_N100_pmax5')
%     merge_hpc_results('S2_N1000_pmax5', true)   % clean up partials
%
%     % Merge all configs at once:
%     configs = {'S1_N100_pmax5', 'S1_N100_pmax10', 'S2_N50_pmax5'};
%     for i = 1:length(configs)
%         merge_hpc_results(configs{i});
%     end

    if nargin < 2, delete_partials = false; end

    % ------------------------------------------------------------------ %
    %  1. Find partial files                                               %
    % ------------------------------------------------------------------ %
    script_dir  = fileparts(mfilename('fullpath'));
    results_dir = fullfile(script_dir, 'results', 'task3');
    partial_dir = fullfile(results_dir, 'partial');
    final_path  = fullfile(results_dir, [config_tag '.mat']);

    pattern     = fullfile(partial_dir, [config_tag '_batch*.mat']);
    batch_files = dir(pattern);

    if isempty(batch_files)
        error('No partial files found for config: %s\nSearched: %s', ...
              config_tag, pattern);
    end

    % Sort by batch number
    [~, sort_idx] = sort({batch_files.name});
    batch_files   = batch_files(sort_idx);
    n_found       = length(batch_files);

    fprintf('[merge] Config: %s\n', config_tag);
    fprintf('[merge] Found %d partial files\n', n_found);

    % ------------------------------------------------------------------ %
    %  2. Load first file to get structure                                 %
    % ------------------------------------------------------------------ %
    d1          = load(fullfile(partial_dir, batch_files(1).name));
    method_lbls = fieldnames(d1.raw);
    n_methods   = length(method_lbls);
    metric_nms  = {'sensitivity','specificity','MCC','FM','HD', ...
                   'TP','FP','TN','FN'};
    gc_true     = d1.gc_true;
    dgp_cfg     = d1.dgp_cfg;
    sim_cfg     = d1.sim_cfg;

    % Determine total MC count
    n_mc_total = 0;
    for f = 1:n_found
        tmp         = load(fullfile(partial_dir, batch_files(f).name), 'batch_info');
        n_mc_total  = n_mc_total + tmp.batch_info.n_mc_batch;
    end
    fprintf('[merge] Total MC runs: %d\n', n_mc_total);

    % ------------------------------------------------------------------ %
    %  3. Initialise merged raw struct                                     %
    % ------------------------------------------------------------------ %
    raw = struct();
    for m = 1:n_methods
        lbl = method_lbls{m};
        for mn = 1:length(metric_nms)
            raw.(lbl).(metric_nms{mn}) = NaN(n_mc_total, 1);
        end
        raw.(lbl).runtime    = NaN(n_mc_total, 1);
        raw.(lbl).n_failures = 0;
    end

    % ------------------------------------------------------------------ %
    %  4. Fill from each partial file in mc_index order                   %
    % ------------------------------------------------------------------ %
    for f = 1:n_found
        fpath  = fullfile(partial_dir, batch_files(f).name);
        d      = load(fpath);
        idxs   = d.batch_info.mc_indices;   % global indices for this batch
        ni     = length(idxs);

        for m = 1:n_methods
            lbl = method_lbls{m};
            if ~isfield(d.raw, lbl), continue; end
            for mn = 1:length(metric_nms)
                if isfield(d.raw.(lbl), metric_nms{mn})
                    raw.(lbl).(metric_nms{mn})(idxs) = ...
                        d.raw.(lbl).(metric_nms{mn})(1:ni);
                end
            end
            raw.(lbl).runtime(idxs)  = d.raw.(lbl).runtime(1:ni);
            raw.(lbl).n_failures     = raw.(lbl).n_failures + d.raw.(lbl).n_failures;
        end
        fprintf('[merge] Loaded batch %d (%d runs)\n', ...
            d.batch_info.batch_id, ni);
    end

    % ------------------------------------------------------------------ %
    %  5. Compute aggregate statistics (same as run_mc_simulation.m)      %
    % ------------------------------------------------------------------ %
    agg = struct();
    for m = 1:n_methods
        lbl = method_lbls{m};
        for mn = 1:length(metric_nms)
            vals = raw.(lbl).(metric_nms{mn});
            vals = vals(~isnan(vals));
            agg.(lbl).(metric_nms{mn}).mean = mean(vals);
            agg.(lbl).(metric_nms{mn}).std  = std(vals);
            agg.(lbl).(metric_nms{mn}).n    = length(vals);
        end
    end

    % ------------------------------------------------------------------ %
    %  6. Print summary table                                             %
    % ------------------------------------------------------------------ %
    fprintf('\n--- Results: %s ---\n', config_tag);
    fprintf('%-12s %6s %6s %6s %6s %6s\n', 'Method','SENS','SPEC','MCC','FM','HD');
    fprintf('%s\n', repmat('-',1,50));
    for m = 1:n_methods
        lbl = method_lbls{m};
        fprintf('%-12s %6.3f %6.3f %6.3f %6.3f %6.2f\n', lbl, ...
            agg.(lbl).sensitivity.mean, agg.(lbl).specificity.mean, ...
            agg.(lbl).MCC.mean, agg.(lbl).FM.mean, agg.(lbl).HD.mean);
    end

    % ------------------------------------------------------------------ %
    %  7. Save final .mat                                                  %
    % ------------------------------------------------------------------ %
    save(final_path, 'raw', 'agg', 'gc_true', 'dgp_cfg', 'sim_cfg', '-v7.3');
    fprintf('\n[merge] Final results saved: %s\n', final_path);

    % ------------------------------------------------------------------ %
    %  8. Optionally delete partial files                                  %
    % ------------------------------------------------------------------ %
    if delete_partials
        for f = 1:n_found
            delete(fullfile(partial_dir, batch_files(f).name));
        end
        fprintf('[merge] Deleted %d partial files.\n', n_found);
    end
end
