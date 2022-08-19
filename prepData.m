%% Here we prepare the model to use.

cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling
load('C:/Work/MatlabCode/components/human-GEM/Human-GEM_1_12/Human-GEM/model/Human-GEM.mat')

%Remove some reactions we do not need - amino acid triplets and drug reactions
m = removeDrugReactions(ihuman);
AATriplets = getAATripletReactions(m,false);
m = removeReactions(m, AATriplets);

%tic
ltModelOrig = CreateECLightModel(m, true, 1);
%toc %Elapsed time is 103.032560 seconds.

ltModelCorr = curateModel(ltModelOrig);
ltModelCorr.lb(ltModelCorr.lb == -1000) = -Inf; %These operations help the solver, it runs faster with inf
ltModelCorr.ub(ltModelCorr.ub == 1000) = Inf;

bloodData = prepBloodData(false, true);
ltModelMin = minimizeModel(ltModelCorr, bloodData); %13790 rxns
save('data/ltModelMin.mat', 'ltModelMin');

cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"

ltModel = setGrowthMedium(ltModelMin, true, 'Hams');
ltModel.lb(strcmp(ltModel.rxns,'MAR03964')) = cell_maintenance;%HMR_3964
save('data/ltModel.mat', 'ltModel');

ltModelFull = setGrowthMedium(ltModelCorr, true, 'Hams');
ltModelFull.lb(strcmp(ltModelFull.rxns,'HMR_3964')) = cell_maintenance;
save('data/ltModelFull.mat', 'ltModelFull');

