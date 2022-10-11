%Generates the data for figure 3 and associated supplementary figures
cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling

%load the models
load('data/ltModel.mat');

cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"

bloodData = prepBloodData();

%Generate a combination of runs with fracECM = 0.01, 0.25, 0.5 - fracC = 0.6, fracF = 0.2
%Then vary frac fibroblast = 0.05, 0.2, 0.4 - fracC is adapted, fracECM = 0.25
%Note that the middle in both are the same, so we get 5 runs in total
%Run each with both type of runs

m1 = addCollaborationCost(buildFullModel2(ltModel, 0.6, 0.2, 0.01, 0.2));
m2 = addCollaborationCost(buildFullModel2(ltModel, 0.6, 0.2, 0.25, 0.2));
m3 = addCollaborationCost(buildFullModel2(ltModel, 0.6, 0.2, 0.5, 0.2));
m4 = addCollaborationCost(buildFullModel2(ltModel, 0.75, 0.05, 0.25, 0.2));
m5 = addCollaborationCost(buildFullModel2(ltModel, 0.45, 0.35, 0.25, 0.2));
m6 = addCollaborationCost(buildFullModel2(ltModel, 0.9699, 0.0001, 0.00001, 0.2));
m6_2 = addCollaborationCost(buildFullModel2(ltModel, 0.7999, 0.0001, 0.00001, 0.2));

%Setup the metabolites for the B runs
%only run for the metabolites of interest since it is heavy to run

selectedExchRxnMets =     {   ...
    'NEFA blood pool in', ...
    'glucose', ...
    'histidine', ...
    'isoleucine', ...
    'leucine', ...
    'lysine', ...
    'methionine', ...
    'phenylalanine', ...
    'threonine', ...
    'tryptophan', ...
    'valine', ...
    'O2', ...
    'alanine', ...
    'asparagine', ...
    'glutamine', ...
    'tyrosine', ...
    'cysteine', ...
    'arginine', ...
    'glycine', ...
    'proline', ...
    'serine', ...
    'aspartate', ...
    'glutamate', ...
    'pyruvate', ...
    'L-lactate', ...
    'albumin'};

%get the metabolites for each exchange reaction
[~, exchRxnInd] = getExchangeRxns(m1, 'in'); %so, the exchange reactions should have the same index in the small and combined model.
exchRxnMets = cell(length(exchRxnInd), 1);

%same here, the mets should be the same
for i = 1:length(exchRxnMets)
   exchRxnMets{i} = m1.metNames{m1.S(:, exchRxnInd(i)) == 1};
end

%now filter so we only get the ones we are interested in.
sel = ismember(exchRxnMets, selectedExchRxnMets);
exchRxnMets = exchRxnMets(sel);
exchRxnInd = exchRxnInd(sel);

%Create mapping between the blood mets and the exch mets:
mappingExchMets = NaN(length(exchRxnMets),1);
for i = 1:length(mappingExchMets)
    ind = find(strcmp(strcat(exchRxnMets(i),'[e]'), bloodData.totMets));
    if length(ind) == 1
        mappingExchMets(i) = ind;
    end
end
%test
tmp = strcat(exchRxnMets(1:(length(exchRxnMets)-1)),'[e]');
all(strcmp(tmp, bloodData.totMets(mappingExchMets(1:(length(exchRxnMets)-1)))))%ok, all equal






%Always use the "standard" range
%a = (0.000001:0.000005:0.0003);
a = (0.000001:0.000001:0.0001); %fewer points since it takes so long to run
aFigA = (0.000001:0.0000015:0.00015); %need longer range for this plot

disp('1A:')
D2_1A = runASimulationFullModel(m1, aFigA, bloodData, cell_maintenance);
save('data/D2_1A.mat', 'D2_1A')

%disp('1B:')
%D2_1B = runMetaboliteImportanceSimulationFullModel(m1, a, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, cell_maintenance, false);
%R2_A = repairASimulationSmallModel(R2_A, ltModel, a, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, cell_maintenance, false);
%save('data/D2_1B.mat', 'D2_1B');

disp('2A:')
D2_2A = runASimulationFullModel(m2, aFigA, bloodData, cell_maintenance);
save('data/D2_2A.mat', 'D2_2A')

disp('3A:')
D2_3A = runASimulationFullModel(m3, aFigA, bloodData, cell_maintenance);
save('data/D2_3A.mat', 'D2_3A')

disp('3B:')
D2_3B = runMetaboliteImportanceSimulation(m3, aFigA, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, cell_maintenance, false);
%R2_A = repairASimulationSmallModel(R2_A, ltModel, a, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, cell_maintenance, false);
save('data/D2_3B.mat', 'D2_3B');


disp('4A:')
D2_4A = runASimulationFullModel(m4, aFigA, bloodData, cell_maintenance);
save('data/D2_4A.mat', 'D2_4A')

disp('5A:')
D2_5A = runASimulationFullModel(m5, aFigA, bloodData, cell_maintenance);
save('data/D2_5A.mat', 'D2_5A')

%Run collaboration experiments (on the m2 model, somewhat in the middle params)
%D2_6tmp = findCollaborationMets(m1, aFigA(100), bloodData, cell_maintenance);
%find(D2_6.collaborationMets)
%for manual test
%inModel = m2
%a = aFigA(60)
%cellMaintenance = cell_maintenance
%and then just run the code in the function line by line

disp('6:')
D2_6 = findCollaborationMets(m1, aFigA, bloodData, cell_maintenance);
save('data/D2_6.mat', 'D2_6');
disp('6B:')
D2_6B = findCollaborationMets(m2, aFigA, bloodData, cell_maintenance);
save('data/D2_6B.mat', 'D2_6B');

disp('7:')
D2_7 = runASimulationFullModel(blockCollaboration(m1), aFigA, bloodData, cell_maintenance);
save('data/D2_7.mat', 'D2_7');

disp('7B:')
D2_7B = runASimulation(ltModel, aFigA, bloodData, cell_maintenance);
save('data/D2_7B.mat', 'D2_7B');

disp('8:')
D2_8 = runASimulationFullModel(m6, aFigA, bloodData, cell_maintenance);
save('data/D2_8.mat', 'D2_8')
%constructEquations(m1, 'o_HMR_9063_REV')
%figure out which substrate main the other cells are using
sel = startsWith(m6.rxns, 'o_');
sum(sel)
hasSomeFlux = D2_8.resultSolutions{20}.x > 10^-8;
sum(hasSomeFlux)
osComp = find(strcmp(m6.comps, 'o_e'));%19
osMets = m6.metComps == osComp;
sComp = find(strcmp(m6.comps, 'e'));%19
sMets = m6.metComps == sComp;
sum(osMets)%1448
exchRxns = (sum(m6.S(osMets,:) ~= 0, 1) > 0) & (sum(m6.S(sMets,:) ~= 0, 1) > 0);
sum(exchRxns)
constructEquations(m6, m6.rxns(hasSomeFlux.' & exchRxns))
%sum them all up across the entire range
%For various values of a, we find that the other cells use different substrates
%for a=20, typically amino acids, citrate, etc.
%for a=50, linoleate

disp('8A:')
D2_8A = runASimulationFullModel(m6_2, aFigA, bloodData, cell_maintenance);
save('data/D2_8A.mat', 'D2_8A')


%% Macrophage simulation
%
load('data/ltModelFull.mat');
disp('9:')

%The macrophage simulation
macrData = readtable('data/MacrophageInput.txt');
D2_9 = runMacrophageSimulation(ltModelFull, a, bloodData, macrData, 0.1); %important to run on the full model - there are new exchange reactions used here that were not used when minimizing the model
save('data/D2_9.mat', 'D2_9')

%And a "normal" run for comparison
disp('10:')
D2_10 = runASimulation(ltModelFull, a, bloodData, cell_maintenance);
save('data/D2_10.mat', 'D2_10')

%% Limit collaboration mets to those mentioned in literature and see what the growth rate is
%
disp('11:')
literatureCollabMets = {'L-lactate','pyruvate','acetone','acetoacetate','(R)-3-hydroxybutanoate','(10Z)-heptadecenoic acid', 'glutamine', 'alanine'};
m1OnlyLiteratureCollab = blockCollaboration(m1, literatureCollabMets);
D2_11 = runASimulationFullModel(m1OnlyLiteratureCollab, aFigA, bloodData, cell_maintenance);
save('data/D2_11.mat', 'D2_11');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the blood flow model
%%%%%%%%%%%%%%%%%%%%%%%%%%%
bloodData2 = prepBloodData(false, false, true);
%test
bloodData2.totDxC(strcmp(bloodData2.totMets,'O2[e]'))

%get the metabolites for each exchange reaction
[~, exchRxnIndAll] = getExchangeRxns(ltModel, 'in');
exchRxnMetsAll = cell(length(exchRxnIndAll), 1);

for i = 1:length(exchRxnMetsAll)
   exchRxnMetsAll{i} = ltModel.metNames{ltModel.S(:, exchRxnIndAll(i)) == 1};
end

%filter mets and exchange rxns that does not exist in the blood data (it is meaningless to work with them)
%filters 5 metabolites + prot_pool
sel = ismember(strcat(exchRxnMetsAll,'[e]'), bloodData2.totMets);
exchRxnMetsAll(~sel)%H2O, Pi, sulfate, Fe2+, hypoxanthine, prot_pool
%check the other way around - all exist there
%sel2 = ismember(bloodData.totMets, strcat(exchRxnMets,'[s]'));
%bloodData.totMets(~sel2)%this is just water etc.



exchRxnMets = exchRxnMetsAll(sel);
exchRxnInd = exchRxnIndAll(sel);

%Create mapping between the blood mets and the exch mets:
mappingExchMets = NaN(length(exchRxnMets),1);
for i = 1:length(mappingExchMets)
    ind = find(strcmp(strcat(exchRxnMets(i),'[e]'), bloodData2.totMets));
    if length(ind) == 1
        mappingExchMets(i) = ind;
    end
end
%test
tmp = strcat(exchRxnMets(1:(length(exchRxnMets)-1)),'[e]');
all(strcmp(tmp, bloodData2.totMets(mappingExchMets(1:(length(exchRxnMets)-1)))))%ok, all equal

load('data/ltModelFull.mat');
disp('9_bloodflow:')




%The "standard" range is roughly 10 times higher here
a_bf = (0.00001:0.00001:0.001);
aFigA_bf = (0.00001:0.000015:0.0015); %need longer range for this plot


%The macrophage simulation
macrData = readtable('data/MacrophageInput.txt');
D2_9_bloodflow = runMacrophageSimulation(ltModelFull, a_bf, bloodData2, macrData, 0.1); %important to run on the full model - there are new exchange reactions used here that were not used when minimizing the model
save('data/D2_9_bloodflow.mat', 'D2_9_bloodflow')

%And a "normal" run for comparison
disp('10_bloodflow:')
D2_10_bloodflow = runASimulation(ltModelFull, a_bf, bloodData2, cell_maintenance);
save('data/D2_10_bloodflow.mat', 'D2_10_bloodflow')

%% Limit collaboration mets to those mentioned in literature and see what the growth rate is
%

%Normal run without collaboration
disp('7_bloodflow:')
D2_7_bloodflow = runASimulationFullModel(blockCollaboration(m1), aFigA_bf, bloodData2, cell_maintenance);
save('data/D2_7_bloodflow.mat', 'D2_7_bloodflow');

%Normal run without collaboration
disp('11:')
literatureCollabMets = {'L-lactate','pyruvate','acetone','acetoacetate','(R)-3-hydroxybutanoate','(10Z)-heptadecenoic acid', 'glutamine', 'alanine'};
m1OnlyLiteratureCollab = blockCollaboration(m1, literatureCollabMets);
D2_11_bloodflow = runASimulationFullModel(m1OnlyLiteratureCollab, aFigA_bf, bloodData2, cell_maintenance);
save('data/D2_11_bloodflow.mat', 'D2_11_bloodflow');


