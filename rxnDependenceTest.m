%Runs the simulations used in for example Fig. 1B, one simulation per value in the a range.
function res = rxnDependenceTest(inModel, aHypox, aNormal, bloodData, cellMaintenance)
%for test
%aHypox = a(30)
%aNormal = a(100)
%inModel = ltModel;
%cellMaintenance = cell_maintenance

% Strategy:
% 1. Run in hypoxia and check which reactions are used. These are the only ones we need to check.
%    Also record the growth rate, both for hypoxic and normal conditions.
% 2. Loop through each reaction and block it, one at the time. Simulate and
%    compare growth rates for hypoxic and normal.

biomassRxn = 'MAR13082';

params = struct();
params.relGap = 0.4;
params.FeasibilityTol = 1e-9;
params.OptimalityTol = 1e-9;

%prepare the hypoxic and normal models
hm = inModel;
ux = zeros(length(bloodData.totDxC),2);
ux(:,2) = bloodData.totDxC*aHypox;
%insert constraints from blood conc + Ham's media
hm = constrainMedium(hm, bloodData.totMets, ux, false, true);
hm.c = double(strcmp(hm.rxns, biomassRxn));%biomass_human set objective function to growth

nm = inModel;
ux = zeros(length(bloodData.totDxC),2);
ux(:,2) = bloodData.totDxC*aNormal;
%insert constraints from blood conc + Ham's media
nm = constrainMedium(nm, bloodData.totMets, ux, false, true);
nm.c = double(strcmp(nm.rxns, biomassRxn));%biomass_human set objective function to growth

%simulate hypoxic and normal
res.baseResH = solveLP(hm,1,params);
res.baseResN = solveLP(nm,1,params);


%Now find the reactions to examine. Only look at reactions with flux in hypoxia and skip 
%exchange reactions:
[~, exchangeRxnsIndexes] = getExchangeRxns(hm);
exchFlt = true(length(hm.rxns), 1);
exchFlt(exchangeRxnsIndexes) = false;
flx = abs(res.baseResH.x) > 0;
rxnSel = exchFlt & flx;
sum(rxnSel) 

res.hyp = cell(length(rxnSel), 1);
res.norm = cell(length(rxnSel), 1);

tot = sum(rxnSel);%771, not so many
curr = 1;
for i = 1:length(rxnSel)
    if rxnSel(i)
        disp([num2str(curr) ' of ' num2str(tot)])
        curr = curr + 1;
        %hypoxic model with blocked rxn
        thm = hm;
        thm.ub(i) = 0;
        thm.lb(i) = 0;
        res.hyp{i} = solveLP(thm,1,params);
        %hypoxic model with blocked rxn
        tnm = nm;
        tnm.ub(i) = 0;
        tnm.lb(i) = 0;
        res.norm{i} = solveLP(tnm,1,params);
    end
end

end