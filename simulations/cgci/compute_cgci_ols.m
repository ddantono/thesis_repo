function [CGCI, pval] = compute_cgci_ols(Y, support, p)
%COMPUTE_CGCI_OLS  Compute OLS-based CGCI for all ordered pairs (i->j).
%
%  For each response j and cause i (i ~= j):
%    1. U-model: OLS regression of Xj on all selected support regressors.
%    2. R-model: same but excluding ALL lags of variable i.
%    3. CGCI(j,i) = log(SSE_R / SSE_U)             -- Eq. 3, Siggiridou 2016
%    4. Fisher statistic (adapted Eq. 8, Siggiridou 2016):
%       F = [(SSE_R - SSE_U)/q_ij] / [SSE_U/df_U]
%       where q_ij = number of i's lags in U-model,
%             df_U  = T_eff - P_j  (T_eff = T-p, P_j = total U-model regressors).
%
%  INPUTS
%    Y        [T x K]     Multivariate time series.
%    support  [K*p x K]  Binary matrix. support(row, col) = 1 means the
%                         regressor at 'row' is included in equation 'col'.
%                         Same convention as the framework B matrix:
%                           row = (lag-1)*K + var_index.
%    p        Positive integer -- VAR lag order. Must satisfy K*p == size(support,1).
%
%  OUTPUTS
%    CGCI  [K x K]  CGCI values.  CGCI(j,i) = CGCI from Xi to Xj.
%                   Diagonal = 0. NaN if df_U <= 0 (over-parameterized).
%    pval  [K x K]  p-values from Fisher test. Diagonal = NaN.
%                   pval = 1 if cause i has no lags in U-model.
%
%  NOTES
%    - No intercept is used (zero-mean VAR assumption, consistent with paper).
%    - Requires Statistics and Machine Learning Toolbox for fcdf().
%    - For dense methods (BGR, HLAG), support = ones(K*p, K) gives
%      standard full-VAR CGCI, matching the "Full" benchmark in the paper.

[T, K] = size(Y);
assert(size(support,1) == K*p && size(support,2) == K, ...
    'compute_cgci_ols: support must be [K*p x K] = [%d x %d].', K*p, K);
assert(T > p, 'compute_cgci_ols: T must be > p.');

T_eff    = T - p;
X_lag    = build_lag_matrix(Y, p);         % [T_eff x K*p]
Y_resp   = Y(p+1:end, :);                  % [T_eff x K]
support  = logical(support);

CGCI = zeros(K, K);
pval = nan(K, K);   % diagonal stays NaN

for j = 1:K   % response equation
    y_j    = Y_resp(:, j);             % [T_eff x 1]
    U_cols = find(support(:, j));      % selected regressor indices (subset of 1..K*p)

    % Fit U-model once per equation (reused for all causes i)
    P_j = length(U_cols);
    df_U = T_eff - P_j;

    if P_j > 0
        X_U = X_lag(:, U_cols);
        if df_U <= 0
            % Over-parameterized — cannot compute CGCI for this equation
            warning('compute_cgci_ols: df_U=%d<=0 for equation j=%d (P_j=%d, T_eff=%d). Setting to NaN.', ...
                    df_U, j, P_j, T_eff);
            for i = 1:K
                if i ~= j
                    CGCI(j,i) = NaN;
                    pval(j,i) = NaN;
                end
            end
            continue;
        end
        b_U   = X_U \ y_j;
        SSE_U = sum((y_j - X_U * b_U).^2);
    else
        % No regressors selected — null model
        SSE_U = sum(y_j.^2);
    end

    for i = 1:K   % cause variable
        if i == j, continue; end

        % All lag indices for variable i in the full regressor matrix
        cols_i = i + (0:p-1)*K;              % [1 x p] indices

        % Lags of i actually present in the U-model
        cols_i_in_U = intersect(U_cols, cols_i);
        q_ij        = length(cols_i_in_U);

        if q_ij == 0
            % Variable i not in U-model -> no evidence of causality
            CGCI(j, i) = 0;
            pval(j, i) = 1;
            continue;
        end

        % R-model: U-model without i's lags
        R_cols = setdiff(U_cols, cols_i_in_U);

        if ~isempty(R_cols)
            X_R   = X_lag(:, R_cols);
            b_R   = X_R \ y_j;
            SSE_R = sum((y_j - X_R * b_R).^2);
        else
            % Null R-model (no remaining regressors)
            SSE_R = sum(y_j.^2);
        end

        % CGCI value (Eq. 3)
        if SSE_U <= 0
            CGCI(j, i) = 0;
            pval(j, i) = 1;
            continue;
        end
        CGCI(j, i) = log(max(SSE_R, SSE_U) / SSE_U);   % clamp: CGCI >= 0

        % Adapted Fisher statistic (Eq. 8)
        F_num = (SSE_R - SSE_U) / q_ij;
        F_den = SSE_U / df_U;

        if F_den <= 0 || F_num <= 0
            pval(j, i) = 1;
        else
            F_stat     = F_num / F_den;
            pval(j, i) = 1 - fcdf(F_stat, q_ij, df_U);
        end
    end
end
end