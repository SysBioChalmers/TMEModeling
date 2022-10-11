%Runs the macrophage simulation, where some of the cells are assumed to die, thereby supplying
%resources for growth.
function res = runMacrophageSimulation(inModel, a, bloodData, macrData, fracReuse)

nPoints = length(a);
res.a = a;
res.resultSolutions = cell(nPoints,1);
res.rxns = inModel.rxns;

for i = 1:nPoints
    disp(i)

    ux = zeros(length(bloodData.totDxC),2);
    ux(:,2) = bloodData.totDxC*a(i);
    %insert constraints from blood conc + Ham's media
    modelGrowth = constrainMedium(inModel, bloodData.totMets, ux, false, true);
    modelGrowth.c = double(strcmp(modelGrowth.rxns,'MAR13082'));%set objective function to growth
    resTmp = solveLP(modelGrowth,1);
    if resTmp.stat ~= -1
        %now add reused materials collected by macrophages
        %first, get the biomass flux from the first simulation
        growth = resTmp.x(strcmp(inModel.rxns,'MAR13082'));%biomass_human
        multFactor = growth * fracReuse; %so, we assume that we reuse for example 10% of the growth, i.e. that those cells die and are collected by the macrophages.
        additionalUb = macrData.TotContent .* multFactor;
        metNames = macrData.Metabolite;
        modelGrowth2 = increaseMetConstraints(modelGrowth, metNames, additionalUb, false, false);
        %test
        %sel = strcmp(modelGrowth2.rxns, 'MAR09033_REV');%'NEFA blood pool in' exch rxn
        %tmp = modelGrowth2.ub(sel) - modelGrowth.ub(sel);
        %expRes = multFactor * 0.5552;
        %expRes - tmp % -2.3419e-17, ok, only a small numerical error

        res.resultSolutions{i} = solveLP(modelGrowth2,1);
    else
        res.resultSolutions{i} = resTmp;
    end
end
