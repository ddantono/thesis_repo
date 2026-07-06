function [raw, agg] = run_mc_simulation(dgp_cfg, method_cfgs, sim_cfg)
%RUN_MC_SIMULATION  Monte Carlo simulation for one (DGP, N, pmax) configuration.
%
%  [raw, agg] = run_mc_simulation(dgp_cfg, method_cfgs, sim_cfg)
%
%  INPUTS
%    dgp_cfg       Struct describing the DGP:
%      .name         char     'S1','S2','S3','S4'
%      .N            int      time series length
%      .pmax         int      VAR lag order for methods and CGCI
%      .K            int      (required for S4) number of variables
%      .C            float    (optional, S4 only) coupling strength (default 0.5)
%      .seed_system  int      (optional, S3 only) seed for B_true topology (default 42)
%
%    method_cfgs   Struct array from define_method_configs().
%
%    sim_cfg       Simulation control:
%      .N_MC         int      number of Monte Carlo runs (default 100)
%      .alpha_fdr    float    FDR level (default 0.05)
%      .burnin       int      burn-in length (default 500)
%      .save_path    char     folder to save results ('' = don't save)
%      .verbose      logical  print per-run progress (default true)
%
%  OUTPUTS
%    raw   Struct: raw.(label).metrics(mc) for each method and MC run.
%    agg   Struct: agg.(label).(metric_name).mean / .std aggregated.
%
%  KEY DESIGN: For S3, B_true is generated ONCE before the loop and
%  reused across all MC runs. Only the innovations vary between runs.

% =========================================================================
%  Defaults
% =========================================================================
if ~isfield(sim_cfg, 'N_MC'),      sim_cfg.N_MC      = 100;   end
if ~isfield(sim_cfg, 'alpha_fdr'), sim_cfg.alpha_fdr = 0.05;  end
if ~isfield(sim_cfg, 'burnin'),    sim_cfg.burnin     = 500;   end
if ~isfield(sim_cfg, 'save_path'), sim_cfg.save_path  = '';    end
if ~isfield(sim_cfg, 'verbose'),   sim_cfg.verbose    = true;  end
if ~isfield(dgp_cfg, 'seed_system'), dgp_cfg.seed_system = 42; end
if ~isfield(dgp_cfg, 'C'),           dgp_cfg.C           = 0.5; end

N_MC     = sim_cfg.N_MC;
N        = dgp_cfg.N;
pmax     = dgp_cfg.pmax;
n_meth   = length(method_cfgs);
metric_names = {'sensitivity','specificity','MCC','FM','HD','TP','FP','TN','FN'};

% =========================================================================
%  Pre-generate fixed DGP topology (critical for S3)
% =========================================================================
fprintf('\n[MC] DGP=%s  N=%d  pmax=%d  N_MC=%d  n_methods=%d\n', ...
        dgp_cfg.name, N, pmax, N_MC, n_meth);

switch upper(dgp_cfg.name)
    case 'S1'
        [~, B_fixed, gc_true] = generate_S1(N, 0, sim_cfg.burnin);
        K = 5;
    case 'S2'
        [~, B_fixed, gc_true] = generate_S2(N, 0, sim_cfg.burnin);
        K = 4;
    case 'S3'
        % Generate B_true ONCE with fixed seed_system
        fprintf('[MC] Generating S3 topology (seed_system=%d)...\n', dgp_cfg.seed_system);
        [~, B_fixed, gc_true] = generate_S3(N, dgp_cfg.seed_system, 0, sim_cfg.burnin);
        K = 20;
    case 'S4'
        K        = dgp_cfg.K;
        B_fixed  = [];   % nonlinear — no B_true
        [~, gc_true] = generate_S4(10, K, dgp_cfg.C, 0, sim_cfg.burnin);
    otherwise
        error('run_mc_simulation: Unknown DGP name: %s', dgp_cfg.name);
end

fprintf('[MC] K=%d,  true GC links=%d / %d\n', K, sum(gc_true(:)), K*(K-1));

% =========================================================================
%  Initialise storage
% =========================================================================
for m = 1:n_meth
    lbl = method_cfgs(m).label;
    raw.(lbl) = struct();
    for mn = 1:length(metric_names)
        raw.(lbl).(metric_names{mn}) = nan(N_MC, 1);
    end
    raw.(lbl).runtime    = nan(N_MC, 1);
    raw.(lbl).n_failures = 0;
end

% =========================================================================
%  Monte Carlo loop
% =========================================================================
t_start = tic;

for mc = 1:N_MC
    % Progress
    if sim_cfg.verbose && (mc == 1 || mod(mc, max(1, floor(N_MC/10))) == 0)
        elapsed = toc(t_start);
        if mc > 1
            eta = elapsed / (mc-1) * (N_MC - mc + 1);
            fprintf('[MC] Run %4d / %d  |  elapsed %.1fs  |  ETA %.1fs\n', ...
                    mc, N_MC, elapsed, eta);
        else
            fprintf('[MC] Run %4d / %d\n', mc, N_MC);
        end
    end

    % --- Generate data for this MC run ---
    seed_noise = mc;   % deterministic seed per run for reproducibility

    switch upper(dgp_cfg.name)
        case {'S1', 'S2'}
            % B_true is fixed; only noise seed varies
            rng(seed_noise, 'twister');
            Y_mc = simulate_var(B_fixed, eye(K), N, sim_cfg.burnin);
        case 'S3'
            % B_true is pre-generated; noise seed varies
            rng(seed_noise, 'twister');
            Y_mc = simulate_var(B_fixed, eye(K), N, sim_cfg.burnin);
        case 'S4'
            [Y_mc, ~] = generate_S4(N, K, dgp_cfg.C, seed_noise, sim_cfg.burnin);
    end

    % --- Run each method ---
    for m = 1:n_meth
        lbl    = method_cfgs(m).label;
        result = run_single_experiment(Y_mc, gc_true, method_cfgs(m), pmax, sim_cfg.alpha_fdr);

        if result.success
            for mn = 1:length(metric_names)
                fname = metric_names{mn};
                raw.(lbl).(fname)(mc) = result.metrics.(fname);
            end
            raw.(lbl).runtime(mc) = result.runtime;
        else
            raw.(lbl).n_failures = raw.(lbl).n_failures + 1;
            if sim_cfg.verbose
                fprintf('  [WARN] MC=%d Method=%s FAILED: %s\n', mc, lbl, result.msg);
            end
        end
    end
end

elapsed_total = toc(t_start);
fprintf('[MC] Complete. Total time: %.1fs  (%.2fs per run)\n', ...
        elapsed_total, elapsed_total/N_MC);

% =========================================================================
%  Aggregate results
% =========================================================================
agg = struct();
for m = 1:n_meth
    lbl = method_cfgs(m).label;
    for mn = 1:length(metric_names)
        fname        = metric_names{mn};
        vals         = raw.(lbl).(fname);
        valid        = ~isnan(vals);
        agg.(lbl).(fname).mean = mean(vals(valid));
        agg.(lbl).(fname).std  = std(vals(valid));
        agg.(lbl).(fname).n    = sum(valid);
    end
    agg.(lbl).runtime.mean   = mean(raw.(lbl).runtime, 'omitnan');
    agg.(lbl).n_failures     = raw.(lbl).n_failures;
end

% =========================================================================
%  Save results
% =========================================================================
if ~isempty(sim_cfg.save_path)
    if ~exist(sim_cfg.save_path, 'dir')
        mkdir(sim_cfg.save_path);
    end
    fname = sprintf('%s_N%d_pmax%d', dgp_cfg.name, N, pmax);
    if strcmp(upper(dgp_cfg.name), 'S4')
        fname = sprintf('%s_K%d_N%d', dgp_cfg.name, K, N);
    end
    save_file = fullfile(sim_cfg.save_path, [fname '.mat']);
    save(save_file, 'raw', 'agg', 'dgp_cfg', 'sim_cfg', 'gc_true');
    fprintf('[MC] Results saved to: %s\n', save_file);
end

% =========================================================================
%  Print summary table
% =========================================================================
print_summary_table(agg, method_cfgs, dgp_cfg);

end


% =========================================================================
%  LOCAL: print results table
% =========================================================================
function print_summary_table(agg, method_cfgs, dgp_cfg)
fprintf('\n--- Results: DGP=%s  N=%d  pmax=%d ---\n', ...
        dgp_cfg.name, dgp_cfg.N, dgp_cfg.pmax);
fprintf('%-12s  %6s  %6s  %6s  %6s  %6s\n', ...
        'Method', 'SENS', 'SPEC', 'MCC', 'FM', 'HD');
fprintf('%s\n', repmat('-', 1, 50));
for m = 1:length(method_cfgs)
    lbl = method_cfgs(m).label;
    a   = agg.(lbl);
    fprintf('%-12s  %6.3f  %6.3f  %6.3f  %6.3f  %6.2f\n', ...
            lbl, ...
            a.sensitivity.mean, a.specificity.mean, ...
            a.MCC.mean, a.FM.mean, a.HD.mean);
end
fprintf('\n');
end