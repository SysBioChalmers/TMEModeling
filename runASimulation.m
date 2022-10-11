%Runs the simulations used in for example Fig. 1B, one simulation per value in the a range.
function res = runASimulation(inModel, a, bloodData, cellMaintenance, relaxOxygen, blockLactateOutputLimit, biomassRxn)
if nargin < 5
    relaxOxygen = false;
end
if nargin < 6
    blockLactateOutputLimit = NaN;
end

if nargin < 7
    biomassRxn = 'MAR13082';
end


nPoints = length(a);
res.a = a;
res.resultSolutions = cell(nPoints,1);
res.rxns = inModel.rxns;
tic

params = struct();
params.relGap = 0.4;
params.FeasibilityTol = 1e-9;
params.OptimalityTol = 1e-9;


for i = 1:nPoints
    disp(i)

    ux = zeros(length(bloodData.totDxC),2);
    ux(:,2) = bloodData.totDxC*a(i);
    %insert constraints from blood conc + Ham's media
    modelGrowth = constrainMedium(inModel, bloodData.totMets, ux, false, true);
    if ~isnan(blockLactateOutputLimit)
        modelGrowth.ub(strcmp(modelGrowth.rxns, 'MAR09135')) = a(i)*blockLactateOutputLimit; % HMR_9135 blockLactateOutputLimit is DxC
    end
    modelGrowth.c = double(strcmp(modelGrowth.rxns, biomassRxn));%biomass_human set objective function to growth
    if relaxOxygen
        modelGrowth.ub(strcmp(modelGrowth.rxns,'MAR09048_REV')) = Inf; %HMR_9048_REV oxygen exchange reaction
    end
    res.resultSolutions{i} = solveLP(modelGrowth,1,params);
end
