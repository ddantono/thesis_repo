function [Y, info] = simulate_var(B, noise_cov, T, burnin)
%SIMULATE_VAR  Simulate a VAR(p) process.
%
%  USAGE
%    [Y, info] = simulate_var(B, noise_cov, T)
%    [Y, info] = simulate_var(B, noise_cov, T, burnin)
%
%  INPUTS
%    B          [K*p x K]  Coefficient matrix (framework convention).
%               B((l-1)*K + i, j) = coeff of variable i at lag l in equation j.
%    noise_cov  Innovation covariance: [K x K] positive-definite matrix,
%               OR scalar sigma^2  (-> sigma^2 * I_K).
%    T          Number of observations to return.
%    burnin     Burn-in length (default: 500).
%
%  OUTPUTS
%    Y          [T x K]  Simulated time series (post burn-in).
%    info       Struct: K, p, T, burnin, spectral_radius, is_stable.
%
%  NOTE
%    The model is: x(t) = A_1*x(t-1) + ... + A_p*x(t-p) + e(t)
%    where A_l(j,i) = B((l-1)*K + i, j)  and  e(t) ~ N(0, noise_cov).

if nargin < 4 || isempty(burnin)
    burnin = 500;
end

% --- Validate B ---
[Kp, K] = size(B);
assert(mod(Kp, K) == 0, ...
    'simulate_var: B must have K*p rows. Got Kp=%d, K=%d.', Kp, K);
p = Kp / K;

% --- Noise covariance ---
if isscalar(noise_cov)
    noise_cov = noise_cov * eye(K);
end
assert(isequal(size(noise_cov), [K K]), ...
    'simulate_var: noise_cov must be [%d x %d] or a scalar.', K, K);

% --- Build A_l matrices ---
% A_l is [K x K]:  A_l(j, i) = B((l-1)*K + i, j)
% Equivalently:    A_l = B(rows_l, :)'  where rows_l = (l-1)*K + (1:K)
A = cell(1, p);
for l = 1:p
    rows_l  = (l-1)*K + (1:K);
    A{l}    = B(rows_l, :)';    % [K x K]
end

% --- Generate innovations ---
T_sim   = T + burnin;
noise   = mvnrnd(zeros(1, K), noise_cov, T_sim);   % [T_sim x K]

% --- Simulate ---
% buf(1:p,:)        initial conditions (near-zero)
% buf(p+1:p+T_sim)  simulated values
buf = zeros(p + T_sim, K);
buf(1:p, :) = randn(p, K) * 1e-3;   % small random start

for t = p+1 : p+T_sim
    xt = zeros(K, 1);
    for l = 1:p
        xt = xt + A{l} * buf(t-l, :)';
    end
    buf(t, :) = xt' + noise(t-p, :);
end

% --- Discard burn-in ---
Y = buf(p + burnin + 1 : p + burnin + T, :);

% --- Build info struct ---
[is_stable, rho] = check_var_stationarity(B);
info = struct('K', K, 'p', p, 'T', T, 'burnin', burnin, ...
              'spectral_radius', rho, 'is_stable', is_stable);

if ~is_stable
    warning(['simulate_var: Process is NOT stationary (rho = %.6f). '  ...
            'Data may diverge.'], rho);
end
end