function metrics = compute_gc_metrics(gc_est, gc_true)
%COMPUTE_GC_METRICS  Evaluate GC network estimation accuracy.
%
%  metrics = compute_gc_metrics(gc_est, gc_true)
%
%  Computes the five performance indices used in Siggiridou & Kugiumtzis (2016),
%  Section II.E, comparing an estimated GC network against the ground truth.
%  The DIAGONAL of both matrices is EXCLUDED from all computations.
%
%  INPUTS
%    gc_est   [K x K]  Estimated binary GC matrix (output of apply_fdr).
%                      gc_est(i,j) = 1 means Xj is estimated to cause Xi.
%    gc_true  [K x K]  True binary GC matrix (from DGP).
%                      gc_true(i,j) = 1 means Xj truly Granger-causes Xi.
%
%  OUTPUT
%    metrics  Struct with fields:
%      TP           True  positives (correctly detected links)
%      FP           False positives (spurious links)
%      TN           True  negatives (correctly absent links)
%      FN           False negatives (missed true links)
%      n_pairs      Total off-diagonal pairs = K*(K-1)
%      n_true       Number of true links
%      sensitivity  TP/(TP+FN)   — recall of true links       [0,1], ideal=1
%      specificity  TN/(TN+FP)   — recall of true non-links   [0,1], ideal=1
%      precision    TP/(TP+FP)   — positive predictive value  [0,1], ideal=1
%      MCC          Matthews Correlation Coefficient           [-1,1], ideal=1
%      FM           F-measure = 2*precision*sensitivity/(precision+sensitivity)
%                                                              [0,1], ideal=1
%      HD           Hamming Distance = FP + FN                [0,K*(K-1)], ideal=0
%
%  EDGE-CASE CONVENTIONS:
%    sensitivity = NaN if TP+FN = 0  (no true links in gc_true)
%    specificity = NaN if TN+FP = 0  (no true non-links in gc_true)
%    precision   = 0   if TP+FP = 0  (nothing detected; consistent with FM=0)
%    MCC         = 0   if denominator = 0  (degenerate confusion matrix)
%    FM          = 0   if precision + sensitivity = 0

K = size(gc_est, 1);
assert(isequal(size(gc_est), size(gc_true)), ...
    'compute_gc_metrics: gc_est and gc_true must be the same size.');
assert(size(gc_est,1) == size(gc_est,2), ...
    'compute_gc_metrics: matrices must be square.');

% Off-diagonal mask (exclude self-pairs)
off_diag = ~logical(eye(K));

% Flatten to vectors (off-diagonal only)
e = gc_est(off_diag);    % estimated, off-diagonal
t = gc_true(off_diag);   % true,      off-diagonal

% Confusion matrix entries
TP = sum( e == 1 &  t == 1);
FP = sum( e == 1 &  t == 0);
TN = sum( e == 0 &  t == 0);
FN = sum( e == 0 &  t == 1);

n_pairs = K * (K - 1);
n_true  = TP + FN;

% ---- Sensitivity (recall) ----
if TP + FN > 0
    sensitivity = TP / (TP + FN);
else
    sensitivity = NaN;   % no true links exist
end

% ---- Specificity ----
if TN + FP > 0
    specificity = TN / (TN + FP);
else
    specificity = NaN;   % no true non-links exist
end

% ---- Precision ----
if TP + FP > 0
    precision = TP / (TP + FP);
else
    precision = 0;   % nothing detected; define precision = 0 for FM continuity
end

% ---- MCC (Matthews Correlation Coefficient) ----
denom_mcc = sqrt(double(TP+FP) * double(TP+FN) * double(TN+FP) * double(TN+FN));
if denom_mcc > 0
    MCC = (double(TP)*double(TN) - double(FP)*double(FN)) / denom_mcc;
else
    MCC = 0;   % degenerate case
end

% ---- F-measure ----
if precision + sensitivity > 0 && ~isnan(sensitivity)
    FM = 2 * precision * sensitivity / (precision + sensitivity);
else
    FM = 0;
end

% ---- Hamming Distance ----
HD = FP + FN;

% ---- Pack output ----
metrics.TP          = TP;
metrics.FP          = FP;
metrics.TN          = TN;
metrics.FN          = FN;
metrics.n_pairs     = n_pairs;
metrics.n_true      = n_true;
metrics.sensitivity = sensitivity;
metrics.specificity = specificity;
metrics.precision   = precision;
metrics.MCC         = MCC;
metrics.FM          = FM;
metrics.HD          = HD;
end