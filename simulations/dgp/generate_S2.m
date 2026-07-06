function [Y, B_true, gc_true] = generate_S2(N, seed, burnin)
%GENERATE_S2  Generate data from DGP S2: sparse VAR(5), K=4.
%
%  Source: Siggiridou & Kugiumtzis (2016) IEEE TSP 64(7).
%  Original model: Winterhalder et al. (2005) Signal Processing 85(11),
%                  2137-2160, model 1.
%
%  SYSTEM EQUATIONS (unit-variance independent innovations ε_i):
%    X1(t) =  0.80·X1(t-1) + 0.65·X2(t-4)                     + ε1
%    X2(t) =  0.60·X2(t-1) + 0.60·X4(t-5)                     + ε2
%    X3(t) =  0.50·X3(t-3) - 0.60·X1(t-1) + 0.40·X2(t-4)     + ε3
%    X4(t) =  1.20·X4(t-1) - 0.70·X4(t-2)                     + ε4
%
%  NOTE ON VAR ORDER:
%    Maximum non-zero lag = 5 (from X4(t-5) in eq.X2). Hence VAR(5).
%    X3's self-coupling is at lag 3 only (lags 1 and 2 are zero for X3).
%    X4 is an exogenous driver: no cross-variable inputs, but it drives X2.
%
%  TRUE CONNECTIVITY (4 directed links):
%    X2->X1  |  X4->X2  |  X1->X3, X2->X3
%    (X4 drives X2 which drives X1 and X3 — a cascade chain)
%
%  SIMULATION PARAMETERS (from paper Table III):
%    Test with pmax = 5 (true order), N = 50, 100, 1000, MC = 1000 runs.
%
%  USAGE
%    [Y, B_true, gc_true] = generate_S2(N)
%    [Y, B_true, gc_true] = generate_S2(N, seed)
%    [Y, B_true, gc_true] = generate_S2(N, seed, burnin)

if nargin < 2 || isempty(seed),   seed   = 0;   end
if nargin < 3 || isempty(burnin), burnin = 500;  end

K = 4;
p = 5;    % VAR(5): maximum non-zero lag is 5 (from X4(t-5) in eq.X2)

% =========================================================================
%  TRUE GC MATRIX  [K x K]
%  gc_true(i,j) = 1  <->  Xj Granger-causes Xi
% =========================================================================
gc_true = zeros(K, K);
gc_true(1, 2) = 1;   % X2 -> X1
gc_true(2, 4) = 1;   % X4 -> X2
gc_true(3, 1) = 1;   % X1 -> X3
gc_true(3, 2) = 1;   % X2 -> X3

% =========================================================================
%  COEFFICIENT MATRIX  [K*p x K] = [20 x 4]
%  Row index: (lag-1)*K + var_index  (1-based, K=4)
%  Column:    equation index
% =========================================================================
B_true = zeros(K*p, K);

% ---- Equation 1 : X1(t) ----
B_true((1-1)*K + 1, 1) =  0.80;  % X1(t-1)
B_true((4-1)*K + 2, 1) =  0.65;  % X2(t-4)   <- X2->X1

% ---- Equation 2 : X2(t) ----
B_true((1-1)*K + 2, 2) =  0.60;  % X2(t-1)
B_true((5-1)*K + 4, 2) =  0.60;  % X4(t-5)   <- X4->X2  [max lag: makes this VAR(5)]

% ---- Equation 3 : X3(t) ----
B_true((3-1)*K + 3, 3) =  0.50;  % X3(t-3)   [no self-lag at 1 or 2]
B_true((1-1)*K + 1, 3) = -0.60;  % X1(t-1)   <- X1->X3
B_true((4-1)*K + 2, 3) =  0.40;  % X2(t-4)   <- X2->X3

% ---- Equation 4 : X4(t) ----  [no cross-variable terms]
B_true((1-1)*K + 4, 4) =  1.20;  % X4(t-1)
B_true((2-1)*K + 4, 4) = -0.70;  % X4(t-2)

% =========================================================================
%  Consistency check: B_true support must match gc_true exactly
% =========================================================================
gc_from_B = zeros(K, K);
for j = 1:K
    rows_j = j + (0:p-1)*K;
    for i = 1:K
        if i ~= j && any(B_true(rows_j, i) ~= 0)
            gc_from_B(i, j) = 1;
        end
    end
end
if ~isequal(gc_from_B, gc_true)
    error('generate_S2: B_true support does not match gc_true. Check equations.');
end

% =========================================================================
%  Stationarity check
% =========================================================================
[is_stable, rho] = check_var_stationarity(B_true);
if ~is_stable
    error('generate_S2: B_true is NOT stationary (rho=%.4f).', rho);
end
fprintf('generate_S2: spectral radius = %.4f  [STABLE]\n', rho);

% =========================================================================
%  Simulate
% =========================================================================
rng(seed, 'twister');
noise_cov = eye(K);
Y = simulate_var(B_true, noise_cov, N, burnin);
end