%Curates the model by blocking some reactions etc.
function [ltModelCorr] = curateModel(ltModelOrig, remFattyAcidRxns)

if nargin < 2
    remFattyAcidRxns = true;
end
%for debugging
%load('C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/ecHumanGEM_batch.mat')
%ecModelOrig = ecModel_batch;

ltModelCorr = ltModelOrig;
%remove some unwanted functions;
%so, we should constrain the reactions HMR_8759, RE1342C_REV and HMR_3996 to have
%zero flux
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR08759')) = 0;%HMR_8759
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR02779_REV')) = 0;%RE1342C_REV
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR03996')) = 0; %HMR_3996 redox path via cystine/cysteine, no GPRs
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR06965')) = 0; %HMR_6965 Reduction path via spermidine. Here there is a GR rule, but many other reactions in this chain are missing GR rules, so this one gets falsely selected regardless
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR13081')) = 0; %CYOOm3i There are two almost identical complex 4 reactions; one with ROS and one without - they are the same reaction. We remove the one with ROS.
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR08981')) = 0; %MDHx Malate to OOA in perixome. No GPR -> bypasses the citric acid cycle reaction, which reduces the proteomic cost of that cycle.
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR03167')) = 0; %RE2591C Production of urate radicals. This function is spontaneous, but is in practice most likely not a main pathway used, so shut it down.
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR12019')) = 0; %LACROX Production of radicals, most likely not the main solution to getting redox power. Shut it down.
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR08759')) = 0; %HMR_8759 This reaction creates a loop where NAD+ can be regenerated without protein cost due to missing GPRs. Shut it down.
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR01706')) = 0; %HMR_1706 This reaction creates a loop in bile metabolism where NADP+ can be regenerated. There is a protein cost, but I suspect this should not be active. Shut it down.
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR08561')) = 0; %HMR_8561 Parallel to p450LTB4r, but not using any NAPDH. looks suspicious, remove.
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR07701')) = 0; %HMR_7701 Regeneration of NADP+ by H2O2 generation, such high levels would probably kill the cell
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR01575')) = 0; %P450LTB4r Redox power regeneration loop from oxygen
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR04423')) = 0; %HMR_4423 Extracellular ROS generation from ornithine=> import and generation of NAD+
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR08606')) = 0; %HMR_8606 Intracellular ROS generation from ornithine=> import and generation of NAD+
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR00059')) = 0; %ACOX22x forms a loop to regenerate NADP+ from oxygen
%ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR03238')) = 0; %RE2675C Regenerates NAD+ from serine and fatty acids (NEFA) - this should probably not be removed
%ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR00754')) = 0; %HMR_0754 Regenerates NADP+ from serine and fatty acids (NEFA) - this should probably not be removed
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR07703')) = 0; %HMR_7703 High ROS production, not realistic
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR07706')) = 0; %HMR_7706 High ROS production, not realistic
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR06539')) = 0; %HMR_6539 High oxygen flux, not realistic
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR06606')) = 0; %HMR_6606 Loop for NAD+ generation including urate radicals, seems unlikely to have a large flux there
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR06611')) = 0; %HMR_6611 Generates H2O2, high flux not realistic (would probably kill the cell)
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR08017')) = 0; %HMR_8017 Generates H2O2, high flux not realistic (would probably kill the cell)
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR08021')) = 0; %HMR_8021 Generates H2O2, high flux not realistic (would probably kill the cell)
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR05390')) = 0; %HMR_5390 Generates H2O2, high flux not realistic (would probably kill the cell)
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR06770')) = 0; %HMR_6770 Generates H2O2, high flux not realistic (would probably kill the cell)
ltModelCorr.ub(strcmp(ltModelCorr.rxns, 'MAR07647')) = 0; %HMR_7647 Generates H2O2, high flux not realistic (would probably kill the cell)


%Temp fix of a bug in the ec model generation - maybe not needed when Gecko
%light is used
%m00205c, m00205g, m00205l, m00205r
ltModelCorr.metNames(strcmp(ltModelCorr.mets, 'MAM00205c')) = {'[protein]-L-serine'};
ltModelCorr.metNames(strcmp(ltModelCorr.mets, 'MAM00205g')) = {'[protein]-L-serine'};
ltModelCorr.metNames(strcmp(ltModelCorr.mets, 'MAM00205l')) = {'[protein]-L-serine'};
ltModelCorr.metNames(strcmp(ltModelCorr.mets, 'MAM00205r')) = {'[protein]-L-serine'};
