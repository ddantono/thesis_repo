function xlagM = multilagmatrix(xM, responseindex, ordersV, indexV)
% MULTILAGMATRIX  Builds the explanatory variable matrix for a DR model.
%
%   xlagM = multilagmatrix(xM, responseindex, ordersV, indexV)
%
%   Called internally by DRfitmse.m and mBTS.m.
%   Do NOT call directly from the framework.
%
%   INPUT
%     xM            : [T x K] matrix of K time series (variables in columns)
%     responseindex : index of response variable in {1,...,K}
%     ordersV       : [1 x K] maximum lag order for each variable
%     indexV        : [1 x K*pmax] binary vector of selected lag terms
%                     e.g. position 2*pmax+3 = 1 means y3(t-3) is selected
%
%   OUTPUT
%     xlagM         : [T_eff x (1+n_selected)] matrix — first column is
%                     the response, remaining columns are selected lags
%
%   Reference: Siggiridou & Kugiumtzis, IEEE TSP, Vol 64(7), 2016

    [n, K] = size(xM);
    pmax   = size(indexV, 2) / K;

    xtempM = NaN(n, K*pmax);
    for iK = 1:K
        xtempM(:, (iK-1)*pmax+1 : ordersV(iK)+(iK-1)*pmax) = ...
            lagmatrix(xM(:,iK), 1:ordersV(iK));
    end

    xlagM = [xM(:, responseindex), xtempM(:, find(indexV == 1))]; %#ok<FNDSB>
    xlagM = xlagM(max(ordersV)+1:end, :);
end