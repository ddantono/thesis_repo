function [Y, B_true, gc_true] = generate_S3(N, seed_system, seed_noise, burnin)
%GENERATE_S3  Generate data from DGP S3.
%
%  S3 is a random sparse VAR(3) process on K=20 variables with approximately
%  10% non-zero coefficients, following Siggiridou & Kugiumtzis (2016),
%  Section III.A (citing Basu & Michailidis 2015 for the generation scheme).
%
%  GENERATION PROCEDURE:
%    1. Initialise B_raw = 0  [K*p x K] = [60 x 20].
%    2. Force AR(1) self-coupling: B_raw((l=1, var=i), eq=i) = 1 for all i.
%    3. Randomly assign the remaining non-zero entries so that the total
%       number of non-zeros equals 10% of K*p*K = 120.
%    4. Scale all non-zero coefficients by a single factor alpha chosen via
%       bisection so that spectral_radius(B_true) < 0.95.
%
%  KEY DESIGN CHOICE:
%    seed_system controls the STRUCTURE of B_true (which links exist).
%    seed_noise  controls the innovation sequence.
%    Use the SAME seed_system across all Monte Carlo runs to keep the
%    network topology fixed; vary seed_noise for different realisations.
%
%  USAGE
%    [Y, B_true, gc_true] = generate_S3(N)
%    [Y, B_true, gc_true] = generate_S3(N, seed_system)
%    [Y, B_true, gc_true] = generate_S3(N, seed_system, seed_noise)
%    [Y, B_true, gc_true] = generate_S3(N, seed_system, seed_noise, burnin)
%
%  INPUTS
%    N            Number of observations.
%    seed_system  RNG seed for generating B_true topology (default: 42).
%    seed_noise   RNG seed for innovation sequence       (default: 0).
%    burnin       Burn-in length                         (default: 500).
%
%  OUTPUTS
%    Y        [N  x 20]   Simulated time series.
%    B_true   [60 x 20]   True coefficient matrix.
%    gc_true  [20 x 20]   Binary GC matrix.
%             gc_true(i,j) = 1 means x_j Granger-causes x_i.

if nargin < 2 || isempty(seed_system), seed_system = 42;  end
if nargin < 3 || isempty(seed_noise),  seed_noise  = 0;   end
if nargin < 4 || isempty(burnin),      burnin      = 500; end

K = 20;
p = 3;

n_total   = K * p * K;               % 1200 total coefficient slots
n_nonzero = round(0.1 * n_total);    % 120 target non-zeros

% =========================================================================
%  Build raw coefficient matrix  (RNG: system topology)
% =========================================================================
rng(seed_system, 'twister');

B_raw = zeros(K*p, K);   % [60 x 20]

% Step 1: forced AR(1) diagonal (lag=1, var=i, eq=i) = 1
forced_lin = zeros(1, K);
for i = 1:K
    row_i           = (1-1)*K + i;          % lag=1, var=i  -> row i
    B_raw(row_i, i) = 1.0;
    % Linear index in column-major storage: (col-1)*nrows + row
    forced_lin(i)   = (i-1)*(K*p) + row_i;
end
n_forced = K;   % 20 AR(1) entries already placed

% Step 2: randomly fill remaining (n_nonzero - n_forced) slots
n_extra  = n_nonzero - n_forced;   % 100
all_lin  = 1 : n_total;
pool     = setdiff(all_lin, forced_lin);
perm     = randperm(numel(pool), n_extra);
B_raw(pool(perm)) = 1.0;

% =========================================================================
%  Scale to stationarity: find largest alpha s.t. rho(alpha*B_raw) < 0.95
% =========================================================================
alpha  = i_find_stability_scale(B_raw, 0.95);
B_true = alpha * B_raw;

[~, rho] = check_var_stationarity(B_true);
n_nz     = sum(B_true(:) ~= 0);
fprintf(['generate_S3 (seed_system=%d): alpha=%.4f | rho=%.4f | '  ...
        'non-zeros=%d / %d  (%.1f%%)\n'], ...
        seed_system, alpha, rho, n_nz, n_total, 100*n_nz/n_total);

% =========================================================================
%  True GC matrix
%  gc_true(i,j) = 1  iff variable j has at least one non-zero entry
%                     in equation i  AND  j ~= i
% =========================================================================
gc_true = zeros(K, K);
for j = 1:K        % candidate cause
    rows_j = j + (0:p-1)*K;    % rows in B for variable j, lags 1..p
    for i = 1:K    % response equation
        if i == j, continue; end
        if any(B_true(rows_j, i) ~= 0)
            gc_true(i, j) = 1;
        end
    end
end

% =========================================================================
%  Simulate  (RNG: noise sequence)
% =========================================================================
rng(seed_noise, 'twister');
noise_cov = eye(K);
Y = simulate_var(B_true, noise_cov, N, burnin);
end


% =========================================================================
%  LOCAL HELPER
% =========================================================================

function alpha = i_find_stability_scale(B, target_rho)
%I_FIND_STABILITY_SCALE  Bisection search for largest alpha so that
%  spectral_radius(alpha * B) < target_rho.

% Fast exit if already stable
[~, rho1] = check_var_stationarity(B);
if rho1 < target_rho
    alpha = 1.0;
    return;
end

% Bisection on [0, 1]
lo = 0.0;
hi = 1.0;
for iter = 1:100
    mid = (lo + hi) / 2;
    [~, rho_mid] = check_var_stationarity(mid * B);
    if rho_mid < target_rho
        lo = mid;
    else
        hi = mid;
    end
    if hi - lo < 1e-8
        break;
    end
end
alpha = lo;
end