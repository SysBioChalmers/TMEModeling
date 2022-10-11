%% Exports met names from the model

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

