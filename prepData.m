%% Here we prepare the model to use.

cd C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/
load('C:/Work/MatlabCode/components/human-GEM/Human-GEM/model/Human-GEM.mat')
%tic
ltModelOrig = CreateECLightModel(ihuman, true, 1);
%toc %Elapsed time is 103.032560 seconds.

ltModelCorr = curateModel(ltModelOrig);
ltModelCorr.lb(ltModelCorr.lb == -1000) = -Inf; %These operations help the solver, it runs faster with inf
ltModelCorr.ub(ltModelCorr.ub == 1000) = Inf;

bloodData = prepBloodData(false, true);
ltModelMin = minimizeModel(ltModelCorr, bloodData);
save('data/ltModelMin.mat', 'ltModelMin');

cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"

ltModel = setGrowthMedium(ltModelMin, true, 'Hams');
ltModel.lb(strcmp(ltModel.rxns,'MAR03964')) = cell_maintenance;%HMR_3964
save('data/ltModel.mat', 'ltModel');

ltModelFull = setGrowthMedium(ltModelCorr, true, 'Hams');
ltModelFull.lb(strcmp(ltModelFull.rxns,'HMR_3964')) = cell_maintenance;
save('data/ltModelFull.mat', 'ltModelFull');

