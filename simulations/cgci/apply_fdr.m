function [gc_binary, threshold] = apply_fdr(pval, alpha)
%APPLY_FDR  Benjamini-Hochberg FDR correction for GC significance testing.
%
%  [gc_binary, threshold] = apply_fdr(pval, alpha)
%
%  Implements the Benjamini-Hochberg (1995) procedure to control the False
%  Discovery Rate across all K*(K-1) ordered variable pairs.
%
%  PROCEDURE:
%    1. Collect the m = K*(K-1) off-diagonal p-values.
%    2. Sort ascending: p_(1) <= p_(2) <= ... <= p_(m).
%    3. Find k* = max{k : p_(k) <= alpha * k / m}.
%    4. Reject (declare significant) all hypotheses with p <= p_(k*).
%
%  Reference: Benjamini Y & Hochberg Y (1995). JRSS-B, 57(1), 289-300.
%
%  INPUTS
%    pval     [K x K]  P-value matrix from compute_cgci_ols.
%                      Diagonal entries are ignored (should be NaN).
%    alpha    Significance level (default: 0.05).
%
%  OUTPUTS
%    gc_binary  [K x K]  Binary GC matrix after FDR correction.
%                        gc_binary(j,i) = 1 means Xi significantly Granger-causes Xj.
%                        Diagonal = 0.
%    threshold  Scalar — the effective p-value threshold applied.
%                        0 if no hypothesis was rejected.

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end

K  = size(pval, 1);
m  = K * (K - 1);    % total number of tests (off-diagonal pairs)

% Extract off-diagonal p-values with their (row, col) indices
rows_idx = zeros(m, 1);
cols_idx = zeros(m, 1);
pvals_vec = zeros(m, 1);

idx = 0;
for j = 1:K
    for i = 1:K
        if i == j, continue; end
        idx = idx + 1;
        rows_idx(idx) = j;
        cols_idx(idx) = i;
        pvals_vec(idx) = pval(j, i);
    end
end

% Handle NaN p-values (over-parameterized equations) — treat as p=1
pvals_vec(isnan(pvals_vec)) = 1;

% Benjamini-Hochberg
[p_sorted, sort_idx] = sort(pvals_vec);
bh_thresholds = alpha * (1:m)' / m;   % alpha*k/m for k=1..m

% Find largest k where p_(k) <= alpha*k/m
significant = (p_sorted <= bh_thresholds);
if any(significant)
    k_star    = find(significant, 1, 'last');
    threshold = p_sorted(k_star);
else
    k_star    = 0;
    threshold = 0;
end

% Build binary GC matrix
gc_binary = zeros(K, K);
if k_star > 0
    rejected_idx = sort_idx(1:k_star);   % original indices of rejected hypotheses
    for r = 1:length(rejected_idx)
        orig = rejected_idx(r);
        gc_binary(rows_idx(orig), cols_idx(orig)) = 1;
    end
end
end