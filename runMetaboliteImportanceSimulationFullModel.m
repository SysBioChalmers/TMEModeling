%Runs the metabolite importance simulation (flux variability analysis).
function res = runMetaboliteImportanceSimulationFullModel(inModel, a, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, cellMaintenance, runHalved)

nPoints = length(a);
res.resultSolutionsBasic = cell(nPoints,1);
res.fluxDivUbs = NaN(nPoints,length(exchRxnMets));
res.metsRedGrowth = NaN(nPoints,length(exchRxnMets));
res.mets = exchRxnMets;
res.a = a;
for i = 1:nPoints
    disp(i)
    ux = zeros(length(bloodData.totDxC),2);
    ux(:,2) = bloodData.totDxC*a(i);
    %insert constraints from blood conc + Ham's media
    modelGrowth = constrainMedium(inModel, bloodData.totMets, ux, false, true);
    modelGrowth.c = double(strcmp(modelGrowth.rxns,'total_biomass'));%set objective function to growth
    modelHalve = modelGrowth;
    %first, maximize for growth
    res.resultSolutionsBasic{i} = solveLP(modelGrowth,0);%don't do parsimonious to save some computation time
    %set min flux of growth to the optimal value (minus some small val to avoid solver issues)
    if (res.resultSolutionsBasic{i}.stat == 1)
        modelGrowth = fixFluxes(modelGrowth, res.resultSolutionsBasic{i}, {'total_biomass'}, false, 10^4);

        for m = 1:length(exchRxnMets)
        %for m = 1:2 %temp
            if isnan(mappingExchMets(m)) %special case for protein pool, seems more difficult to solve, use more margin
                modelGrowth = fixFluxes(modelGrowth, res.resultSolutionsBasic{i}, {'total_biomass'}, false, 10^3);
            end
            disp([num2str(i) ':' num2str(m)]);
            %change objective to minimize the metabolite
            modelGrowth.c(:) = 0;
            modelGrowth.c(exchRxnInd(m)) = -1;
            resTmp = solveLP(modelGrowth,0);%don't do parsimonious to save some computation time
            if resTmp.stat == 1
                res.fluxDivUbs(i,m) = resTmp.x(exchRxnInd(m))./modelGrowth.ub(exchRxnInd(m));
            else
                res.fluxDivUbs(i,m) = NaN;
            end

            if (runHalved)
                %Now reduce the metabolite uptake to 90%
                if ~isnan(mappingExchMets(m))
                    uxTmp = ux;
                    uxTmp(mappingExchMets(m),2) = uxTmp(mappingExchMets(m),2)*0.9;
                    modelTmp = constrainMedium(modelHalve, bloodData.totMets, uxTmp, false, true);

                    resTmp = solveLP(modelTmp,0);%don't do parsimonious to save some computation time
                    if resTmp.stat == 1
                        res.metsRedGrowth(i,m) = -resTmp.f;        
                    else
                        res.metsRedGrowth(i,m) = NaN;        
                    end
                end
            end
        end
    end
end
