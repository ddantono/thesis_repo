function [Y, gc_true] = generate_S4(N, K, C, seed, burnin)
%GENERATE_S4  Simulate coupled Hénon maps (DGP S4).
%
%  The map for variable i is:
%
%    x_{i,t} = 1.4 - [ 0.5*C*(x_{i-1,t-1} + x_{i+1,t-1})
%                     + (1-C)*x_{i,t-1} ]^2  +  0.3*x_{i,t-2}
%
%  BOUNDARY CONDITIONS (open chain, missing neighbour treated as zero):
%    i = 1  : no left  neighbour -> x_{0,t-1} = 0
%    i = K  : no right neighbour -> x_{K+1,t-1} = 0
%
%  COUPLING TOPOLOGY (symmetric tridiagonal chain):
%    gc_true(i, i-1) = 1  for  i = 2..K   (left  neighbour drives i)
%    gc_true(i, i+1) = 1  for  i = 1..K-1 (right neighbour drives i)
%
%  Reference: Politi & Torcini (1992); Siggiridou & Kugiumtzis (2016) S4.
%
%  USAGE
%    [Y, gc_true] = generate_S4(N, K)
%    [Y, gc_true] = generate_S4(N, K, C)         % default C = 0.5
%    [Y, gc_true] = generate_S4(N, K, C, seed)
%    [Y, gc_true] = generate_S4(N, K, C, seed, burnin)
%
%  INPUTS
%    N       Number of observations to return.
%    K       Number of variables. Typical values: 5, 10, 20.
%    C       Coupling strength in [0, 1] (default: 0.5).
%    seed    RNG seed for initial conditions (default: 0).
%    burnin  Burn-in length (default: 2000 — longer needed for chaos).
%
%  OUTPUTS
%    Y        [N x K]  Simulated time series (post burn-in).
%    gc_true  [K x K]  Binary GC matrix.
%             gc_true(i,j) = 1  means x_j Granger-causes x_i.
%             Equivalently: gc_true(i,j) = 1  iff  abs(i-j) == 1.
%
%  NOTES
%    - S4 is a NONLINEAR system; no B_true is defined.
%    - Default C=0.5 replicates the coupling used in Siggiridou (2016).
%    - The system is deterministic and chaotic; the only source of
%      variability across realisations is the initial condition (seed).
%    - Use a longer burnin for larger K to ensure the transient decays.

if nargin < 3 || isempty(C),      C      = 0.5;  end
if nargin < 4 || isempty(seed),   seed   = 0;    end
if nargin < 5 || isempty(burnin), burnin = 2000; end

assert(K >= 3,           'generate_S4: K must be >= 3.');
assert(C >= 0 && C <= 1, 'generate_S4: C must be in [0, 1].');
assert(N > 0,            'generate_S4: N must be positive.');

% =========================================================================
%  True GC matrix  (symmetric tridiagonal, diagonal = 0)
%  gc_true(i,j) = 1  iff  |i - j| == 1
% =========================================================================
idx      = (1:K)';
gc_true  = double(abs(bsxfun(@minus, idx, idx')) == 1);

% =========================================================================
%  Simulate
% =========================================================================
rng(seed, 'twister');

T_total = N + burnin;

% Allocate: rows are time, columns are variables
% Need 2 extra rows at the start for lag-2 initialisation.
X = zeros(T_total + 2, K);
X(1:2, :) = rand(2, K) * 0.2 - 0.1;   % small random initial state

for t = 3 : T_total + 2
    for i = 1:K
        % Neighbour values at t-1 (zero for out-of-chain)
        x_left  = 0;
        x_right = 0;
        if i > 1, x_left  = X(t-1, i-1); end
        if i < K, x_right = X(t-1, i+1); end

        % Coupling term
        coupling = 0.5*C*(x_left + x_right) + (1 - C)*X(t-1, i);
        X(t, i)  = 1.4 - coupling^2 + 0.3*X(t-2, i);
    end
end

% Discard burn-in (first 2 rows are initial conditions, not counted)
Y = X(2 + burnin + 1 : 2 + burnin + N, :);

% Sanity check for divergence (can occur with C near boundaries)
if any(~isfinite(Y(:)))
    warning(['generate_S4: Non-finite values detected. '  ...
            'Try a different seed or reduce C.']);
end
end