function res = runASimulationFullModel(inModel, a, bloodData, cellMaintenance, relaxOxygen)
if nargin < 5
    relaxOxygen = false;
end
nPoints = length(a);
res.a = a;
res.resultSolutions = cell(nPoints,1);
res.rxns = inModel.rxns;
tic

for i = 1:nPoints
    disp(i)

    ux = zeros(length(bloodData.totDxC),2);
    ux(:,2) = bloodData.totDxC*a(i);
    %insert constraints from blood conc + Ham's media
    modelGrowth = constrainMedium(inModel, bloodData.totMets, ux, false, true);
    modelGrowth.c = double(strcmp(modelGrowth.rxns,'total_biomass'));%set objective function to cancer cell growth, including ECM
    if relaxOxygen
        modelGrowth.ub(strcmp(modelGrowth.rxns,'MAR09048_REV')) = Inf; %oxygen exchange reaction
    end
    res.resultSolutions{i} = solveLP(modelGrowth,1);
end
