%Helper function to show fluxes from the result of a simulation
function res = listExchRxnsWithFluxes(model, sol, uptake, lowerLimit)
if nargin < 3
    lowerLimit = 0;
end

if uptake
    [~,exchRxnInd] = getExchangeRxns(model, 'in');
else
    [~,exchRxnInd] = getExchangeRxns(model, 'out');
end

fluxSel = sol.x > lowerLimit;
totSel = false(length(sol.x),1);
totSel(exchRxnInd) = true;
totSel = totSel & fluxSel;

fracFluxOfUb = sol.x./model.ub;

tbl = table(model.rxns(totSel), sol.x(totSel), fracFluxOfUb(totSel), constructEquations(model, model.rxns(totSel)));
tbl.Properties.VariableNames = {'Rxn', 'Flux', 'FracFluxOfUb', 'Formula'};
tbl
res = tbl;

end
