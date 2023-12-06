cd C:/Work/MatlabCode/projects/HumanGemZenodo/Human1_Publication_Data_Scripts/ec_GEMs/ComplementaryScripts

cellLines = {'HOP62'};
i=1;

load(['../models/' cellLines{i} '/ecModel_batch.mat'])
load(['../models/' cellLines{i} '/ltModel.mat'])

commonMets = intersect(ecModel_batch.mets, ltModel.mets);
commonRxns = intersect(ecModel_batch.rxns, ltModel.rxns);
%remove prot_pool
commonMets = commonMets(~strcmp(commonMets,'prot_pool'));
ltS = ltModel.S(ismember(ltModel.mets, commonMets), ismember(ltModel.rxns, commonRxns));
ltRxns = ltModel.rxns(ismember(ltModel.rxns, commonRxns));
ltMets = ltModel.mets(ismember(ltModel.mets, commonMets));

ecS = ecModel_batch.S(ismember(ecModel_batch.mets, commonMets), ismember(ecModel_batch.rxns, commonRxns));
ecRxns = ecModel_batch.rxns(ismember(ecModel_batch.rxns, commonRxns));
ecMets = ecModel_batch.mets(ismember(ecModel_batch.mets, commonMets));

all(strcmp(ltRxns,ecRxns)) %ok, in the right order
all(strcmp(ltMets,ecMets)) %ok

%find reactions that are different
diffRxns = false(length(ecRxns),1);
for i = 1:length(ecRxns)
    diffRxns(i,1) = ~all(ltS(:,i) == ecS(:,i));
end
diffRxns
sum(diffRxns) %1942

% ind = find(diffRxns)
% rxn = ltRxns(ind(1))
% constructEquations(ltModel, rxn)
% constructEquations(ecModel_batch, rxn)

% constructEquations(ltModel, 'biomass_human')
% constructEquations(ecModel_batch, 'biomass_human')


ecModel_batch.rxns(find(ecModel_batch.c))
ecModel_batch.ub(strcmp(ecModel_batch.rxns, 'prot_pool_exchange'))
ltModel.ub(strcmp(ltModel.rxns, 'prot_pool_exchange'))

ec = ecModel_batch;
lt = ltModel;

%table(ec.rxns(1:5000),constructEquations(ec,ec.rxns(1:5000)))

%constructEquations(ec,{'HMR_3944No1'})%{'glucose-1-phosphate[c] + H+[c] + UTP[c] + 2.9728e-06 prot_Q16851[s] => PPi[c] + UDP-glucose[c]'}
%constructEquations(ec,{'arm_HMR_3944'})%{'glucose-1-phosphate[c] + H+[c] + UTP[c] + 2.9728e-06 prot_Q16851[s] => PPi[c] + UDP-glucose[c]'}

%{'HMR_4473No1'               }    {'6-phospho-D-gluconate[r] + NADP+[r] + 3.8052e-06 prot_P52209[s] => CO2[r] + NADPH[r] + ribulose-5-phosphate[r]'}
%constructEquations(lt,{'HMR_4473'})%{'glucose-1-phosphate[c] + H+[c] + UTP[c] + 2.9728e-06 prot_Q16851[s] => PPi[c] + UDP-glucose[c]'}
% %{'6-phospho-D-gluconate[r] + NADP+[r] + 0.0002022 prot_pool[c] => CO2[r] + NADPH[r] + ribulose-5-phosphate[r]'}
% met = find(strcmp(ec.metNames,'prot_P52209'));
% met
% rxns = find(ec.S(met,:) ~= 0);
% rxns 
% constructEquations(ec,ec.rxns(rxns)) %{'53.1395 prot_pool[s] => prot_P52209[s]' 
% %is this the same?
% 53.1395 * 3.8052e-06 %2.0221e-04, the same as in GeckoLight

ec2 = setHamsMedium_GeckoLight(ec,1);
lt2 = setHamsMedium_GeckoLight(lt,1);

rec = solveLP(ec2,1);
rlt = solveLP(lt2,1);
rlt

ec2.ub(strcmp(ec2.rxns,'prot_pool_exchange'))
lt2.ub(strcmp(lt2.rxns,'prot_pool_exchange'))

ec2.ub(ec2.ub ~= Inf & ec2.ub ~= 1000 & ec2.ub ~= 0)

s = rec.x ~= 0;
table(ec2.rxns(s), rec.x(s))

constructEquations()






constructEquations(ec,{'HMR_10062'})

HMR_10062
HMR_10063
HMR_10064
HMR_10065
biomass_human








