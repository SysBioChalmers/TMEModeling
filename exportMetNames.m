%% Here we test the full tumor model.

%load the precalculated EC model
%load('C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/ecHumanGEM_batch.mat')
%ecModelOrig = ecModel_batch;

%pmet = startsWith(ecModelOrig.metNames, 'pmet_')
%sum(pmet)%3791
%prot = startsWith(ecModelOrig.metNames, 'prot_')
%sum(prot)%3225

%ecModelOrig.metNames(~(pmet|prot))
cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling
load('C:/Work/MatlabCode/components/human-GEM/Human-GEM/model/Human-GEM.mat')

hirrev = convertToIrrev(ihuman);

exchRxns = getExchangeRxns(hirrev, 'in');
length(exchRxns)%1665
indices = getIndexes(hirrev, exchRxns, 'rxns');
SFilt = hirrev.S(:,indices);
metsSel = sum(SFilt,2) > 0;
hirrev.comps %so 1 == s

metsImportedByModel = hirrev.metNames(metsSel);
save('data/metsImportedByModel.mat', 'metsImportedByModel');

