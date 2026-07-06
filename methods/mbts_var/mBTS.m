function [indexV, maxorder, MSEval] = mBTS(xM, responseindex, pmax)
% MBTS  Modified Backward-in-Time Selection for sparse VAR estimation.
%
%   [indexV, maxorder, MSEval] = mBTS(xM, responseindex, pmax)
%
%   Called equation-by-equation from run_mbts_var.m.
%   Selects only the lag terms found significantly explanatory for the
%   response variable, using BIC as selection criterion.
%
%   INPUT
%     xM            : [T x K] matrix of K time series (variables in columns)
%     responseindex : index of response variable in {1,...,K}
%     pmax          : maximum lag order to search
%
%   OUTPUT
%     indexV        : [1 x K*pmax] binary vector of selected lag terms
%                     in lag-major order (all vars at lag 1, then lag 2, etc.)
%     maxorder      : maximum lag selected across all variables
%     MSEval        : mean squared error of the fitted DR model
%
%   Copyright (C) 2021 Dimitris Kugiumtzis — GPL v3
%   Reference: Siggiridou & Kugiumtzis, IEEE TSP, Vol 64(7), 2016

    [n, K] = size(xM);

    % Remove mean of each time series
    for d = 1:K
        xM(:,d) = xM(:,d) - mean(xM(:,d));
    end

    % Initially no lagged variable is selected
    indexV    = zeros(1, K*pmax);
    ordersinV = zeros(1, K);

    % Initial MSE = variance of response, BIC gets largest value
    MSEval  = DRfitmse(xM, responseindex, ordersinV, indexV);
    BICold  = (n - max(ordersinV)) * log(MSEval) + ...
               sum(ordersinV) * log(n - max(ordersinV));

    ingameV       = 1:K;
    ningame       = K;
    terminateflag = 0;
    incrisor      = 1;

    while (~terminateflag && ningame ~= 0)
        pmaxreach = find(ordersinV >= pmax);
        ingameV   = setdiff(ingameV, pmaxreach);
        ningame   = length(ingameV);

        if ningame ~= 0
            BICnowV = NaN(ningame, 1);
            MSEnowV = NaN(ningame, 1);

            for iK = 1:ningame
                ordtempV = ordersinV;
                ordtempV(ingameV(iK)) = ordtempV(ingameV(iK)) + incrisor;

                overpmaxV = find(ordtempV > pmax);
                if ~isempty(overpmaxV)
                    ordtempV(overpmaxV) = pmax;
                end
                if length(overpmaxV) == K
                    terminateflag = 1;
                end

                tempindexV = indexV;
                tempindexV((ingameV(iK)-1)*pmax + ordtempV(ingameV(iK))) = 1;

                MSEnowV(iK) = DRfitmse(xM, responseindex, ordtempV, tempindexV);
                BICnowV(iK) = (n - max(ordtempV)) * log(MSEnowV(iK)) + ...
                               sum(tempindexV) * log(n - max(ordtempV));
            end

            [BICnew, iBICnew] = min(BICnowV);
            invarindex = ingameV(BICnowV == BICnew);

            if BICold <= BICnew
                incrisor = incrisor + 1;
                if incrisor > pmax - min(ordersinV)
                    terminateflag = 1;
                end
            else
                indexV((invarindex-1)*pmax + ordersinV(invarindex) + incrisor) = 1;
                ordersinV(invarindex) = ordersinV(invarindex) + incrisor;
                BICold   = BICnew;
                incrisor = 1;
                MSEval   = MSEnowV(iBICnew);
            end
        else
            terminateflag = 1;
        end
    end

    % Reshape indexV to lag-major order:
    % all variables at lag 1, then lag 2, etc.
    indexV   = reshape(indexV, pmax, K);
    indexV   = indexV';
    indexV   = reshape(indexV, 1, K*pmax);
    maxorder = max(ordersinV);
end