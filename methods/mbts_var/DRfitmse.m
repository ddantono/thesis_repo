function MSE = DRfitmse(xM, responseindex, ordersV, indexV)
% DRFITMSE  Fits a dynamic regression (DR) model and returns the MSE.
%
%   MSE = DRfitmse(xM, responseindex, ordersV, indexV)
%
%   Called internally by mBTS.m during lag selection.
%   Do NOT call directly from the framework.
%
%   INPUT
%     xM            : [T x K] matrix of K time series (variables in columns)
%     responseindex : index of response variable in {1,...,K}
%     ordersV       : [1 x K] maximum lag order for each variable
%     indexV        : [1 x K*pmax] binary vector of selected lag terms
%                     e.g. position 2*K+3 = 1 means y2(t-3) is selected
%
%   OUTPUT
%     MSE           : mean squared error of the DR model fit
%
%   Reference: Siggiridou & Kugiumtzis, IEEE TSP, Vol 64(7), 2016

    n         = size(xM, 1);
    xlagM     = multilagmatrix(xM, responseindex, ordersV, indexV);
    An        = inv(xlagM(:,2:end)' * xlagM(:,2:end)) * ...
                    xlagM(:,2:end)' * xlagM(:,1);
    preV      = xlagM(:,2:end) * An;
    resV      = xM(max(ordersV)+1:end, responseindex) - preV;
    MSE       = resV' * resV / length(resV);
end