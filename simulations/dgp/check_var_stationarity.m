function [is_stable, rho] = check_var_stationarity(B)
%CHECK_VAR_STATIONARITY  Assess VAR stability via the companion matrix.
%
%  A VAR(p) process is (covariance) stationary iff all eigenvalues of the
%  companion matrix lie strictly inside the unit circle (spectral radius < 1).
%
%  USAGE
%    [is_stable, rho] = check_var_stationarity(B)
%
%  INPUT
%    B  [K*p x K]  VAR coefficient matrix (framework convention).
%
%  OUTPUTS
%    is_stable  logical  — true iff spectral radius < 1
%    rho        scalar   — maximum absolute eigenvalue of the companion matrix
%
%  Companion matrix structure (size K*p x K*p):
%
%    C = [ A_1   A_2   ...  A_p  ]   <- top K rows, each A_l is [K x K]
%        [ I_K   0     ...  0    ]
%        [ 0     I_K   ...  0    ]
%        [ ...                   ]
%        [ 0     0     ...  0    ]
%
%  where  A_l(j,i) = B((l-1)*K + i, j).

[Kp, K] = size(B);
assert(mod(Kp, K) == 0, ...
    'check_var_stationarity: B must have K*p rows. Got Kp=%d, K=%d.', Kp, K);
p = Kp / K;

C = zeros(K*p, K*p);

% Top block: A_l in column-block position l
for l = 1:p
    rows_l      = (l-1)*K + (1:K);
    A_l         = B(rows_l, :)';            % [K x K]
    C(1:K, rows_l) = A_l;
end

% Sub-diagonal identity blocks
if p > 1
    C(K+1 : K*p, 1 : K*(p-1)) = eye(K*(p-1));
end

ev        = eig(C);
rho       = max(abs(ev));
is_stable = (rho < 1 - 1e-10);
end