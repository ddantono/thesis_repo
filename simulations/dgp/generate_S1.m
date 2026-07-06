function [Y, B_true, gc_true] = generate_S1(N, seed, burnin)
%GENERATE_S1  Generate data from DGP S1: sparse VAR(4), K=5.
%
%  Source: Siggiridou & Kugiumtzis (2016) IEEE TSP 64(7), equation (9).
%  Original model: Schelter et al. (2006) J Physiol-Paris 99(1), 37-46.
%
%  ALL EQUATIONS CONFIRMED from paper screenshot.
%
%  SYSTEM EQUATIONS (unit-variance independent innovations u_i):
%
%    X1(t) =  0.4·X1(t-1) - 0.5·X1(t-2) + 0.4·X5(t-1)          + u1
%    X2(t) =  0.4·X2(t-1) - 0.3·X1(t-4) + 0.4·X5(t-2)          + u2
%    X3(t) =  0.5·X3(t-1) - 0.7·X3(t-2) - 0.3·X5(t-3)          + u3
%    X4(t) =  0.8·X4(t-3) + 0.4·X1(t-2) + 0.3·X2(t-2)          + u4
%    X5(t) =  0.7·X5(t-1) - 0.5·X5(t-2) - 0.4·X4(t-1)          + u5
%
%  NOTE ON VAR ORDER:
%    Maximum non-zero lag = 4 (from X1(t-4) in equation X2). VAR(4).
%    X4 has self-coupling only at lag 3 (lags 1 and 2 are zero for X4).
%
%  TRUE CONNECTIVITY (7 directed links, verified from paper p.1766):
%    X5->X1 | X1->X2, X5->X2 | X5->X3 | X1->X4, X2->X4 | X4->X5
%
%  SIMULATION PARAMETERS matching paper:
%    pmax = 5 (slightly over-specified) or pmax = 10 (Table II).
%    N = 100 (short), MC = 1000 runs.
%
%  USAGE
%    [Y, B_true, gc_true] = generate_S1(N)
%    [Y, B_true, gc_true] = generate_S1(N, seed)
%    [Y, B_true, gc_true] = generate_S1(N, seed, burnin)

if nargin < 2 || isempty(seed),   seed   = 0;   end
if nargin < 3 || isempty(burnin), burnin = 500;  end

K = 5;
p = 4;    % VAR(4): maximum non-zero lag = 4 (X1 at lag 4 in eq.X2)

% =========================================================================
%  TRUE GC MATRIX  [K x K]
%  gc_true(i,j) = 1  <->  Xj Granger-causes Xi
% =========================================================================
gc_true = zeros(K, K);
gc_true(1, 5) = 1;   % X5 -> X1
gc_true(2, 1) = 1;   % X1 -> X2
gc_true(2, 5) = 1;   % X5 -> X2
gc_true(3, 5) = 1;   % X5 -> X3
gc_true(4, 1) = 1;   % X1 -> X4
gc_true(4, 2) = 1;   % X2 -> X4
gc_true(5, 4) = 1;   % X4 -> X5

% =========================================================================
%  COEFFICIENT MATRIX  [K*p x K] = [20 x 5]
%
%  Row index formula:  row = (lag - 1)*K + var_index    (1-based, K=5)
%
%    lag 1: rows  1.. 5   (vars 1..5 at lag 1)
%    lag 2: rows  6..10   (vars 1..5 at lag 2)
%    lag 3: rows 11..15   (vars 1..5 at lag 3)
%    lag 4: rows 16..20   (vars 1..5 at lag 4)
%
%  Column = equation index (1..K)
% =========================================================================
B_true = zeros(K*p, K);

% ---- Equation 1 : X1(t) ----
B_true((1-1)*K + 1, 1) =  0.4;   % X1(t-1)  → row  1
B_true((2-1)*K + 1, 1) = -0.5;   % X1(t-2)  → row  6
B_true((1-1)*K + 5, 1) =  0.4;   % X5(t-1)  → row  5   [X5->X1]

% ---- Equation 2 : X2(t) ----
B_true((1-1)*K + 2, 2) =  0.4;   % X2(t-1)  → row  2
B_true((4-1)*K + 1, 2) = -0.3;   % X1(t-4)  → row 16   [X1->X2, max lag → VAR(4)]
B_true((2-1)*K + 5, 2) =  0.4;   % X5(t-2)  → row 10   [X5->X2]

% ---- Equation 3 : X3(t) ----
B_true((1-1)*K + 3, 3) =  0.5;   % X3(t-1)  → row  3
B_true((2-1)*K + 3, 3) = -0.7;   % X3(t-2)  → row  8
B_true((3-1)*K + 5, 3) = -0.3;   % X5(t-3)  → row 15   [X5->X3]

% ---- Equation 4 : X4(t) ----
B_true((3-1)*K + 4, 4) =  0.8;   % X4(t-3)  → row 14   [self at lag 3 only]
B_true((2-1)*K + 1, 4) =  0.4;   % X1(t-2)  → row  6   [X1->X4]
B_true((2-1)*K + 2, 4) =  0.3;   % X2(t-2)  → row  7   [X2->X4]

% ---- Equation 5 : X5(t) ----
B_true((1-1)*K + 5, 5) =  0.7;   % X5(t-1)  → row  5
B_true((2-1)*K + 5, 5) = -0.5;   % X5(t-2)  → row 10
B_true((1-1)*K + 4, 5) = -0.4;   % X4(t-1)  → row  4   [X4->X5]

% =========================================================================
%  Consistency check: B_true non-zero support must match gc_true exactly
% =========================================================================
gc_from_B = zeros(K, K);
for j = 1:K
    rows_j = j + (0:p-1)*K;   % all row indices for variable j
    for i = 1:K
        if i ~= j && any(B_true(rows_j, i) ~= 0)
            gc_from_B(i, j) = 1;
        end
    end
end
if ~isequal(gc_from_B, gc_true)
    error('generate_S1: B_true support does not match gc_true. Check equations.');
end

% =========================================================================
%  Stationarity check
%  Expected: spectral radius dominated by X4 AR(lag-3) component ≈ 0.93
% =========================================================================
[is_stable, rho] = check_var_stationarity(B_true);
if ~is_stable
    error('generate_S1: B_true is NOT stationary (rho=%.4f).', rho);
end
fprintf('generate_S1: spectral radius = %.4f  [STABLE]\n', rho);

% =========================================================================
%  Simulate
% =========================================================================
rng(seed, 'twister');
noise_cov = eye(K);
Y = simulate_var(B_true, noise_cov, N, burnin);
end