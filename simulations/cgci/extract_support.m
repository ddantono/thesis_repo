function support = extract_support(out, tol)
%EXTRACT_SUPPORT  Extract binary support matrix from a framework canonical output.
%
%  support = extract_support(out)
%  support = extract_support(out, tol)
%
%  Returns a binary [K*p x K] matrix where support(row, col) = 1 means the
%  regressor at 'row' is included in equation 'col'.
%
%  Works uniformly across ALL method types:
%    Sparse methods  (LASSO, mBTS, Adaptive LASSO, SCAD, MCP, gglasso):
%      B has exact or near-exact zeros -> support = abs(B) > tol
%
%    Dense methods   (BGR/Ridge, HLAG):
%      B is fully non-zero -> support = all ones
%      OLS-CGCI on full support = standard full-VAR CGCI (the "Full" benchmark)
%
%    Methods with native gc_matrix (pds_lm_var, msvar):
%      B is still in out.coefficients -> same treatment.
%      Their native gc_matrix in out.diagnostics.gc_matrix is a SECONDARY output
%      and can be compared separately; the PRIMARY evaluation uses OLS-CGCI here.
%
%  INPUTS
%    out  Canonical framework output struct (must have field 'coefficients').
%    tol  Threshold below which an entry is treated as zero (default: 1e-10).
%
%  OUTPUT
%    support  [K*p x K]  Logical binary support matrix.

if nargin < 2 || isempty(tol)
    tol = 1e-10;
end

assert(isfield(out, 'coefficients'), ...
    'extract_support: out must have a ''coefficients'' field.');
assert(isnumeric(out.coefficients) && ~isempty(out.coefficients), ...
    'extract_support: out.coefficients must be a non-empty numeric matrix.');

support = abs(out.coefficients) > tol;
end