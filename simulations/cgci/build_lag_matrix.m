function X_lag = build_lag_matrix(Y, p)
%BUILD_LAG_MATRIX  Construct the lagged regressor matrix for a VAR(p) model.
%
%  X_lag = build_lag_matrix(Y, p)
%
%  INPUTS
%    Y      [T x K]  Multivariate time series.
%    p      Positive integer — VAR lag order.
%
%  OUTPUT
%    X_lag  [T_eff x K*p]  Lagged regressor matrix, T_eff = T - p.
%
%  COLUMN CONVENTION (matches framework B matrix row convention):
%    Column (l-1)*K + k  =  variable k at lag l
%    This matches: B_row (l-1)*K + k  =  variable k at lag l
%
%  EXAMPLE
%    For K=2, p=3:
%      col 1 = X1(t-1), col 2 = X2(t-1)
%      col 3 = X1(t-2), col 4 = X2(t-2)
%      col 5 = X1(t-3), col 6 = X2(t-3)

[T, K] = size(Y);
assert(p >= 1, 'build_lag_matrix: p must be >= 1.');
assert(T > p,  'build_lag_matrix: T must be > p.');

T_eff = T - p;
X_lag = zeros(T_eff, K*p);

for l = 1:p
    cols     = (l-1)*K + (1:K);   % columns in X_lag for lag l
    rows_Y   = (p - l + 1) : (T - l);   % rows of Y shifted by lag l
    X_lag(:, cols) = Y(rows_Y, :);
end
end