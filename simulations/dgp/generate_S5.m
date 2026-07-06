function [Y, B_true, gc_true] = generate_S5(N, seed_system, seed_noise, burnin)
%GENERATE_S5  Generate data from DGP S5: sparse VAR(2), K=50.
%
%  S5 is a random sparse VAR(2) process on K=50 variables with approximately
%  10% non-zero coefficients. It is the extreme high-dimensional DGP of
%  this benchmark study, extending the dimensionality of S3 (K=20) by a
%  factor of 2.5.
%
%  MOTIVATION (why S5 uses a different strategy than S3):
%    S3 initialises all non-zero support entries to a uniform value of 1.0
%    and finds a single global scaling factor alpha via bisection. At K=20
%    this yields alpha ≈ 0.15, which is small but still produces detectable
%    causal effects. At K=50 the same strategy forces alpha into the range
%    0.03–0.06, making virtually all couplings indistinguishable from noise
%    and defeating the purpose of the experiment.
%
%    S5 therefore adopts the VARPstable strategy (Kugiumtzis, personal
%    communication): non-zero coefficients are initialised with random
%    heterogeneous values (mixed signs, varied magnitudes) and the entire
%    coefficient matrix is iteratively deflated by 5% per step until the
%    spectral radius of the companion matrix falls below 0.95.  This
%    preserves the relative structure of the coefficient magnitudes during
%    deflation, ensuring that the largest couplings remain meaningfully
%    above the noise floor even in the high-dimensional regime.
%
%  GENERATION PROCEDURE:
%    1. Using seed_system, build a binary support mask on B [K*p x K]:
%       - Force AR(1) self-coupling at lag 1 for all K variables (50 entries).
%       - Randomly fill the remaining 450 non-zero slots from the pool of
%         non-forced positions, targeting 10% density of K*p*K = 5000 slots.
%    2. Using seed_system (still active), initialise each non-zero entry
%       with a heterogeneous random value:
%         sign(0.5 + randn) * (0.3 + 0.9 * rand)
%       The bias 0.5 in the sign term skews toward positive values (~62%),
%       consistent with VARPstable.  Magnitudes range from 0.3 to 1.2.
%    3. Iteratively deflate: B_true = B_true * 0.95 until rho < 0.95.
%       This is an outer loop that calls check_var_stationarity at every
%       iteration.  Convergence is guaranteed because the spectral radius
%       scales continuously with the coefficient magnitudes.
%    4. Derive gc_true from the non-zero support of B_true (off-diagonal).
%    5. Using seed_noise, simulate N observations via simulate_var.
%
%  KEY DESIGN CHOICES:
%    seed_system  Controls the network TOPOLOGY (which links exist) AND the
%                 coefficient magnitudes/signs.  Fix this across all MC runs
%                 to keep the true network constant.
%    seed_noise   Controls the innovation sequence.  Vary this across Monte
%                 Carlo realisations (use the global MC index as seed_noise).
%    burnin       Set to 500 (same as S3) — adequate for a VAR(2) process.
%
%  SYSTEM PARAMETERS:
%    K = 50   (variables)
%    p = 2    (true VAR order; maximum non-zero lag = 2)
%    Total parameter slots : K*p*K = 5000
%    Target non-zeros      : 500  (10%)
%    Forced AR(1) entries  : 50   (one per variable, diagonal at lag 1)
%    Random extra entries  : 450
%    Target spectral radius: rho < 0.95
%
%  USAGE
%    [Y, B_true, gc_true] = generate_S5(N)
%    [Y, B_true, gc_true] = generate_S5(N, seed_system)
%    [Y, B_true, gc_true] = generate_S5(N, seed_system, seed_noise)
%    [Y, B_true, gc_true] = generate_S5(N, seed_system, seed_noise, burnin)
%
%  INPUTS
%    N            Number of observations to return.
%    seed_system  RNG seed for topology and coefficient values (default: 99).
%    seed_noise   RNG seed for innovation sequence            (default: 0).
%    burnin       Burn-in length                              (default: 500).
%
%  OUTPUTS
%    Y        [N  x 50]   Simulated time series (post burn-in).
%    B_true   [100 x 50]  True VAR(2) coefficient matrix [K*p x K].
%                         B_true((l-1)*K + i, j) = coefficient of X_i at
%                         lag l in equation j  (framework convention).
%    gc_true  [50 x 50]   Binary Granger causality matrix.
%                         gc_true(i,j) = 1  means X_j Granger-causes X_i.
%
%  DEPENDENCIES
%    check_var_stationarity.m  (simulations/dgp/)
%    simulate_var.m            (simulations/dgp/)
%
%  REFERENCE
%    Generation strategy based on VARPstable.m (D. Kugiumtzis, personal
%    communication), adapted to the [K*p x K] framework convention and
%    extended with explicit dual-seed control and stationarity targeting.

% =========================================================================
%  Default arguments
% =========================================================================
if nargin < 2 || isempty(seed_system), seed_system = 99;  end
if nargin < 3 || isempty(seed_noise),  seed_noise  = 0;   end
if nargin < 4 || isempty(burnin),      burnin      = 500; end

% =========================================================================
%  Fixed system parameters
% =========================================================================
K          = 50;
p          = 2;
n_total    = K * p * K;           % 5000 total coefficient slots
n_nonzero  = round(0.1 * n_total);% 500  target non-zeros (10% density)
n_forced   = K;                   % 50   forced AR(1) diagonal entries
n_extra    = n_nonzero - n_forced;% 450  additional random non-zero entries

% =========================================================================
%  STAGE 1: Build support and initialise coefficients  (RNG: seed_system)
% =========================================================================
rng(seed_system, 'twister');

B_true = zeros(K*p, K);   % [100 x 50]

% --- Step 1a: Force AR(1) self-coupling for all K variables ---
%  Row index for variable i at lag l: (l-1)*K + i  (1-based, framework conv.)
%  For lag l=1, variable i: row = i.  Column = i (equation i).
forced_lin = zeros(1, K);
for i = 1:K
    row_i              = i;                       % lag=1, var=i -> row i
    % Heterogeneous initialisation for the forced AR(1) entries.
    % sign(0.5 + randn): biased-positive (~62% positive).
    % magnitude: 0.3 + 0.9*rand -> range [0.30, 1.20].
    B_true(row_i, i)   = sign(0.5 + randn()) * (0.3 + 0.9 * rand());
    % Linear index in column-major storage: (col-1)*nrows + row
    forced_lin(i)      = (i - 1) * (K*p) + row_i;
end

% --- Step 1b: Randomly select 450 additional non-zero positions ---
all_lin  = 1 : n_total;
pool     = setdiff(all_lin, forced_lin);
perm     = randperm(numel(pool), n_extra);
extra_positions = pool(perm);

% --- Step 1c: Assign heterogeneous random values to extra positions ---
%  Same distribution as Step 1a: biased-positive signs, varied magnitudes.
for k = 1:n_extra
    B_true(extra_positions(k)) = sign(0.5 + randn()) * (0.3 + 0.9 * rand());
end

% =========================================================================
%  STAGE 2: Iterative deflation to stationarity (VARPstable strategy)
%
%  Multiply the entire B_true by 0.95 at each iteration until
%  spectral_radius(B_true) < 0.95.
%
%  This preserves the relative structure of coefficient magnitudes (unlike
%  the bisection approach which imposes a single global scale on a uniform
%  support), ensuring that strong couplings remain detectably above the
%  noise floor after convergence.
% =========================================================================
[is_stable, rho] = check_var_stationarity(B_true);

n_iter = 0;
while ~is_stable || rho >= 0.95
    B_true   = B_true * 0.95;
    n_iter   = n_iter + 1;
    [is_stable, rho] = check_var_stationarity(B_true);

    % Safety valve: abort after 2000 iterations (should never be needed).
    if n_iter > 2000
        error('generate_S5: deflation did not converge after 2000 iterations (rho=%.4f).', rho);
    end
end

% =========================================================================
%  STAGE 3: Diagnostic report
% =========================================================================
n_nz     = sum(B_true(:) ~= 0);
% After deflation all support entries are still non-zero (deflation
% scales but never zeroes), so n_nz equals n_nonzero by construction.
fprintf(['generate_S5 (seed_system=%d): ' ...
         'deflation_iters=%d | rho=%.4f | ' ...
         'non-zeros=%d / %d (%.1f%%)\n'], ...
         seed_system, n_iter, rho, n_nz, n_total, 100 * n_nz / n_total);

% =========================================================================
%  STAGE 4: Derive true Granger causality matrix
%
%  gc_true(i,j) = 1  iff variable j has at least one non-zero coefficient
%                     in equation i  AND  j ~= i  (no self-links reported).
%
%  Row indices for variable j across all lags:
%    rows_j = j, j+K, j+2K, ..., j+(p-1)*K   (1-based)
%  i.e., rows_j = j + (0:p-1)*K
% =========================================================================
gc_true = zeros(K, K);
for j = 1:K        % candidate driving variable
    rows_j = j + (0:p-1)*K;    % row indices of variable j at lags 1..p
    for i = 1:K    % response equation
        if i == j, continue; end
        if any(B_true(rows_j, i) ~= 0)
            gc_true(i, j) = 1;
        end
    end
end

n_links  = sum(gc_true(:));
n_possible = K * (K - 1);   % 2450 directed non-self pairs
fprintf('generate_S5 (seed_system=%d): GC links = %d / %d (%.1f%%)\n', ...
        seed_system, n_links, n_possible, 100 * n_links / n_possible);

% =========================================================================
%  STAGE 5: Consistency check
%  Verify that the non-zero off-diagonal support of B_true matches gc_true
%  exactly.  This catches any indexing error in the gc_true construction.
% =========================================================================
gc_from_B = zeros(K, K);
for j = 1:K
    rows_j = j + (0:p-1)*K;
    for i = 1:K
        if i == j, continue; end
        if any(B_true(rows_j, i) ~= 0)
            gc_from_B(i, j) = 1;
        end
    end
end
if ~isequal(gc_from_B, gc_true)
    error('generate_S5: internal consistency check failed — B_true support does not match gc_true.');
end

% =========================================================================
%  STAGE 6: Simulate  (RNG: seed_noise)
%  Switch RNG to seed_noise so that the noise sequence is independent of
%  the system topology.  This is the only source of variation across
%  Monte Carlo realisations (B_true and gc_true are fixed by seed_system).
% =========================================================================
rng(seed_noise, 'twister');
noise_cov = eye(K);
Y = simulate_var(B_true, noise_cov, N, burnin);
end