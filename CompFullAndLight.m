%Here, we use the full and light implementations in GECKO 3 to compare full and light versions

adapterLocation = fullfile(findGECKOroot,'tutorials','light_ecModel','HumanGEMAdapter.m');
adapter = ModelAdapterManager.setDefault(adapterLocation); 

model = load('C:/Work/MatlabCode/components/human-GEM/Human-GEM_1_12/Human-GEM/model/Human-GEM.mat').ihuman;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Light
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ecModelLt, noUniprot] = makeEcModel(model,true); %warning is ok
ecModelLt         = getECfromGEM(ecModelLt);
kcatList_fuzzy_lt  = fuzzyKcatMatching(ecModelLt);

ecModelLt = selectKcatValue(ecModelLt,kcatList_fuzzy_lt);
ecModelLt = applyKcatConstraints(ecModelLt);
ecModelLt = setProtPoolSize(ecModelLt);
find(isnan(ecModelLt.S(:,strcmp(ecModelLt.rxns,'MAR08263'))))
%ecModelLt.S(isnan(ecModelLt.S)) = 0;
slt = solveLP(ecModelLt,1)%-0.1428

ecModelLt %8370 mets, 17538 rxns

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Full
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ecModelFull, noUniprot] = makeEcModel(model,false); %warning is ok
ecModelFull         = getECfromGEM(ecModelFull);
kcatList_fuzzy_full  = fuzzyKcatMatching(ecModelFull);

ecModelFull = selectKcatValue(ecModelFull,kcatList_fuzzy_full);
ecModelFull = applyKcatConstraints(ecModelFull);
ecModelFull = setProtPoolSize(ecModelFull);
%ecModelFull.S(isnan(ecModelFull.S)) = 0;
sfull = solveLP(ecModelFull,1)%-0.1443

ecModelFull %11254 mets, 43448 rxns

