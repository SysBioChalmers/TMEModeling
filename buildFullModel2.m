function [outModel] = buildFullModel2(ecModelCorr, fracC, fracF, fracECM, fracGAGsInECM)
% buildFullModel
%
% Generates a combined ec model with tumor cells, fibroblasts, and other cells (immune cells etc).
% The model supports a total tumor biomass objective function which builds
% ECM and tumor cells. The other cell types are assumed to be imported into
% the tumor.
%
% Input:
%
%   ecModelCorr     Should be an ec model derived from Human-GEM
%
%   fracC           fracCion of cancer cells in the model
%
%   fracF           fracCion of fibroblasts in the model. fracCion other cells 
%                   is calculated as 1 - fracC - fracF
%
%   fracECM         fracCion extracellular matrix in the tumor. The rest is
%                   cells
%
%   fracGAGsInECM   fracCion of GAGs in the ECM. The rest is protein.
%
% Output:
%
%   outModel        The full model
%

%for debugging
%load('C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/ecHumanGEM_batch.mat')
%ecModelOrig = ecModel_batch;

%length(getExchangeRxns(ecModelCorr)) %3347
%constructEquations(ecModelCorr, getExchangeRxns(ecModelCorr))
%sum(getMetsInComp(ecModelCorr, 'S'))
%ecModelCorr.metNames(getMetsInComp(ecModelCorr, 'S'))

%ecModelCorr.metComps

%figure out which metabolite is the protein pool, i.e. in the
%'prot_pool_exchange' reaction
%constructEquations(ecModelCorr, ecModelCorr.rxns(length(ecModelCorr.rxns))) %{' => prot_pool[c]'}
%So, this is not an exchange flux from the s compartment, which means we do not need to have any
%special treatment of this when copying the model

%The strategy for creating the three cell types combined model is to
%1. create 3 copies of all metabolites, where the metabolites in the s
%compartment are not copied.
%2. create 3 copies of all reactions, where the exchange reactions are
%excluded
%3. In the matrix, the second copy of the reactions should have the second
%copy of the metabolites. The first and third should be zero, except for
%the [s] metabolites in the first part.


%Identify metabolites to copy
metsToCopy = true(length(ecModelCorr.mets),1);
%identify rxns to copy
rxnsToCopy = true(1,length(ecModelCorr.rxns));

model2 = ecModelCorr;

%generate new compartments
nms = ecModelCorr.comps;
nmsF = strcat('f_', nms);
nmsO = strcat('o_', nms);
model2.comps = [ecModelCorr.comps;nmsF;nmsO];
model2.compNames = [ecModelCorr.compNames;strcat('f_',ecModelCorr.compNames);strcat('o_',ecModelCorr.compNames)];

%extend the metabolites
nms = ecModelCorr.mets(metsToCopy);
nmsF = strcat('f_', nms);
nmsO = strcat('o_', nms);
model2.mets = [ecModelCorr.mets; nmsF; nmsO];
model2.metNames = [ecModelCorr.metNames; ecModelCorr.metNames(metsToCopy); ecModelCorr.metNames(metsToCopy)];
%met compartments
metComps = ecModelCorr.metComps(metsToCopy);

toAdd = length(ecModelCorr.comps);
compsF = metComps + toAdd;
compsO = metComps + toAdd.*2;
%unique(compsF)
%unique(compsO) %looks good
model2.metComps = [ecModelCorr.metComps;compsF;compsO];
%b
model2.b = [ecModelCorr.b;ecModelCorr.b(metsToCopy);ecModelCorr.b(metsToCopy)];

%now, the reactions
nms = ecModelCorr.rxns(rxnsToCopy);
nmsF = strcat('f_', nms);
nmsO = strcat('o_', nms);
model2.rxns = [ecModelCorr.rxns; nmsF; nmsO];
nms = ecModelCorr.rxnNames(rxnsToCopy);
nmsF = strcat('f_', nms);
nmsO = strcat('o_', nms);
model2.rxnNames = [ecModelCorr.rxnNames; nmsF; nmsO];

%lb, ub, rev and c
model2.lb = [ecModelCorr.lb;ecModelCorr.lb(rxnsToCopy);ecModelCorr.lb(rxnsToCopy)];
model2.ub = [ecModelCorr.ub;ecModelCorr.ub(rxnsToCopy);ecModelCorr.ub(rxnsToCopy)];
model2.rev = [ecModelCorr.rev;ecModelCorr.rev(rxnsToCopy);ecModelCorr.rev(rxnsToCopy)];
model2.c = [ecModelCorr.c;ecModelCorr.c(rxnsToCopy);ecModelCorr.c(rxnsToCopy)];

%grRules subSystems and eccodes
model2.grRules = [ecModelCorr.grRules;ecModelCorr.grRules(rxnsToCopy);ecModelCorr.grRules(rxnsToCopy)];
model2.subSystems = [ecModelCorr.subSystems;ecModelCorr.subSystems(rxnsToCopy);ecModelCorr.subSystems(rxnsToCopy)];
model2.eccodes = [ecModelCorr.eccodes;ecModelCorr.eccodes(rxnsToCopy);ecModelCorr.eccodes(rxnsToCopy)];

%rxnGeneMat
model2.rxnGeneMat = [ecModelCorr.rxnGeneMat;ecModelCorr.rxnGeneMat(rxnsToCopy,:);ecModelCorr.rxnGeneMat(rxnsToCopy,:)];

%And now finally the S matrix:
%first just add zeros in both directions
model2.S(length(model2.mets), length(model2.rxns)) = 0;
[nm,nr] = size(ecModelCorr.S);
model2.S(1:nm,1:nr) = ecModelCorr.S;
model2.S((1:nm)+nm,(1:nr)+nr) = ecModelCorr.S;
model2.S((1:nm)+2*nm,(1:nr)+2*nr) = ecModelCorr.S;

%Now, redirect the exchange reactions. The uptake reactions of both f and o
%should take up everything from the S compartment, not from nothing. They
%should also be unconstrained
rxnSelC = 1:nr;
rxnSelF = (1:nr) + nr;
rxnSelO = (1:nr) + nr*2;
metSelC = (1:nm).';
metSelF = (1:nm).' + nm;
metSelO = (1:nm).' + nm*2;

[inExchRxns,inExchRxnsInd] = getExchangeRxns(ecModelCorr, 'in');
%skip protein pool exchange
inExchRxnsInd = inExchRxnsInd(~strcmp(inExchRxns,'prot_pool_exchange'));
%constructEquations(model2,inExchRxns)%looks good
%constructEquations(model2,model2.rxns(inExchRxnsInd + nr))%looks good
model2.S(metSelC, rxnSelF(inExchRxnsInd)) = -ecModelCorr.S(metSelC, rxnSelC(inExchRxnsInd));
model2.S(metSelC, rxnSelO(inExchRxnsInd)) = -ecModelCorr.S(metSelC, rxnSelC(inExchRxnsInd));
%constructEquations(model2,model2.rxns(inExchRxnsInd + nr))%looks good
%constructEquations(model2,model2.rxns(inExchRxnsInd + 2*nr))%looks good

%relax the constraints - the input to the S compartment is constrained
model2.ub(inExchRxnsInd + nr) = Inf;
model2.ub(inExchRxnsInd + 2*nr) = Inf;

%Now redirect the output of the fibroblasts S compartment to the S
%compartment instead of to nothing
[outExchRxns,outExchRxnsInd] = getExchangeRxns(ecModelCorr, 'out');
%constructEquations(ecModelCorr,ecModelCorr.rxns(outExchRxnsInd))%
%constructEquations(model2,model2.rxns(outExchRxnsInd + nr))%
model2.S(metSelC, rxnSelF(outExchRxnsInd)) = -ecModelCorr.S(metSelC, rxnSelC(outExchRxnsInd));
%constructEquations(model2,model2.rxns(outExchRxnsInd + nr))%looks good


%% Now add the specialized objective functions and maintenance

%So, we define a few variables:
%cell type fracCions
%fracC = 0.6; - input param to function
%fracF = 0.2; - input param to function
fracO = 1 - fracC - fracF;
%ECM vs cells fracCion (dry weight)
%fracECM = 0.3; - input param to function
%fracGAGsInECM = 0.2; - input param to function
% We normalize everything around the cancer cells, meaning that we don't
% need to modify the biomass function or the protein pool constraint for
% cancer cells
% This means that the protein pool for fibroblasts should be scaled with 
% fracF/fracC and for other cells fracO/fracC
% It also means that the maintenance cost (i.e. lower bound) should be
% scaled with the same constant.
% The ECM function should be scaled to produce the same dryweight biomass as 
% the tumor cell biomass function scaled with 
% fracECM/((1-fracECM)*fracC)

%Add maintenance
cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"
model2.lb(strcmp(model2.rxns,'MAR03964')) = cell_maintenance; %tumor cells
model2.lb(strcmp(model2.rxns,'f_MAR03964')) = cell_maintenance * fracF/fracC; %fibroblasts
model2.lb(strcmp(model2.rxns,'o_MAR03964')) = cell_maintenance * fracO/fracC; %tumor cells

%change protein pool
t_prot_ub = model2.ub(strcmp(model2.rxns,'prot_pool_exchange'));
model2.ub(strcmp(model2.rxns,'f_prot_pool_exchange')) = t_prot_ub * fracF/fracC; %fibroblasts
model2.ub(strcmp(model2.rxns,'o_prot_pool_exchange')) = t_prot_ub * fracO/fracC; %other cells

%Add reactions for ecm

modelBackup = model2;
model2 = modelBackup;

%first ECM_protein_pool_biomass, i.e. the collagen biomass
metsToAdd.mets = {'ECM_protein_pool_biomass', 'ECM_protein_pool_biomass_e'};
metsToAdd.metNames = {'ECM_protein_pool_biomass', 'ECM_protein_pool_biomass'};
metsToAdd.compartments = {'f_c', 'f_e'};
model2 = addMets(model2, metsToAdd, false);

rxnsToAdd.rxns = {'ecm_protein_pool_biomass'};
rxnsToAdd.equations = {'0.0948 L-alanyl-tRNA(ala)[f_c] + 0.0499 L-arginyl-tRNA(arg)[f_c] + 0.0228 L-asparaginyl-tRNA(asn)[f_c] + 0.0405 L-aspartyl-tRNA(asp)[f_c] + 0.0104 L-cysteinyl-tRNA(cys)[f_c] + 0.0304 L-glutaminyl-tRNA(gln)[f_c] + 0.0503 L-glutamyl-tRNA(glu)[f_c] + 0.2710 glycyl-tRNA(gly)[f_c] + 0.0078 L-histidyl-tRNA(his)[f_c] + 0.0187 L-isoleucyl-tRNA(ile)[f_c] + 0.0367 L-leucyl-tRNA(leu)[f_c] + 0.0382 L-lysyl-tRNA(lys)[f_c] + 0.0084 L-methionyl-tRNA(met)[f_c] + 0.0177 L-phenylalanyl-tRNA(phe)[f_c] + 0.1832 L-prolyl-tRNA(pro)[f_c] + 0.0400 L-seryl-tRNA(ser)[f_c] + 0.0307 L-threonyl-tRNA(thr)[f_c] + 0.0040 L-tryptophanyl-tRNA(trp)[f_c] + 0.0098 L-tyrosyl-tRNA(tyr)[f_c] + 0.0348 L-valyl-tRNA(val)[f_c] => ECM_protein_pool_biomass[f_c] + 0.0948 tRNA(ala)[f_c] + 0.0499 tRNA(arg)[f_c] + 0.0228 tRNA(asn)[f_c] + 0.0405 tRNA(asp)[f_c] + 0.0104 tRNA(cys)[f_c] + 0.0304 tRNA(gln)[f_c] + 0.0503 tRNA(glu)[f_c] + 0.2710 tRNA(gly)[f_c] + 0.0078 tRNA(his)[f_c] + 0.0187 tRNA(ile)[f_c] + 0.0367 tRNA(leu)[f_c] + 0.0382 tRNA(lys)[f_c] + 0.0084 tRNA(met)[f_c] + 0.0177 tRNA(phe)[f_c] + 0.1832 tRNA(pro)[f_c] + 0.0400 tRNA(ser)[f_c] + 0.0307 tRNA(thr)[f_c] + 0.0040 tRNA(trp)[f_c] + 0.0098 tRNA(tyr)[f_c] + 0.0348 tRNA(val)[f_c]'};
model2 = addRxns(model2,rxnsToAdd, 3);
%and, export it out
rxnsToAdd.rxns = {'ecm_protein_pool_biomass_transport'};
rxnsToAdd.equations = {'ECM_protein_pool_biomass[f_c] => ECM_protein_pool_biomass[f_e]'};
model2 = addRxns(model2,rxnsToAdd, 3);

%Now the GAGs
%The strategy is to not include the cost for generating the protein
%structures, but only the xylose and attached heparan sulfate.
%The metabolite to produce is heparan sulfate proteoglycan[s]
%We should shut down possible production of this in other cell types, i.e.
%set the flux (ub) from [g] and [o_g] to [e] to zero in t (HMR_7223). Also shut down 
%degradation of it in the lysosome (HMR_7224).
%To make this work, we also need to provide protein to attach the GAGs on.
%We don't want to create the proteins, since that is accounted for in the 
%ECM_protein_pool_biomass. Instead, we open up an unlimited influx of such
%protein scaffolds, i.e. an exchange reaction of [protein]-L-serine[f_c]
%Also, shut down generation of this protein scaffold in the fibroblasts, to make
%sure it is not used: f_HMR_7197 should have zero flux

%first, shut down heparan sulfate proteoglycan[s] production from other
%cells:
model2.ub(strcmp(model2.rxns,'MAR07223')) = 0;
model2.ub(strcmp(model2.rxns,'o_MAR07223')) = 0;

%and not degraded at all
model2.ub(strcmp(model2.rxns,'MAR07224')) = 0;
model2.ub(strcmp(model2.rxns,'f_MAR07224')) = 0;
model2.ub(strcmp(model2.rxns,'o_MAR07224')) = 0;

%Now, create an exchange reaction for uptake of [protein]-L-serine[f_c]

model2.ub(strcmp(model2.rxns,'f_MAR07197')) = 0; %Make sure protein scaffold is not created inside the cell

%add the protein scaffold in the [e] compartment
metsToAdd.mets = {'gag_scaffold_protein_e'};
metsToAdd.metNames = {'[protein]-L-serine'};
metsToAdd.compartments = {'f_e'};
model2 = addMets(model2, metsToAdd, false);

%and the import reactions
rxnsToAdd.rxns = {'gag_scaffold_protein_import', 'gag_scaffold_protein_exchange'};
rxnsToAdd.equations = {'[protein]-L-serine[f_e] => [protein]-L-serine[f_c]', '=> [protein]-L-serine[f_e]'};
model2 = addRxns(model2,rxnsToAdd, 3);

%Now create an ECM biomass reaction containing both GAGs and protein
protPoolBiomassMW = 94.81;
gagBiomassMW = 3133.44;
protPoolStoichiometry = (1-fracGAGsInECM)/protPoolBiomassMW*1000;
gagStoichiometry = fracGAGsInECM/gagBiomassMW*1000;

metsToAdd.mets = {'ECM_biomass'};
metsToAdd.metNames = {'ECM_biomass'};
metsToAdd.compartments = {'f_e'};
model2 = addMets(model2, metsToAdd, false);

rxn = [num2str(protPoolStoichiometry), ' ECM_protein_pool_biomass[f_e] + ', ...
       num2str(gagStoichiometry), ' heparan sulfate proteoglycan[f_e]', ...
       ' => ECM_biomass[f_e]'];
%rxn = [num2str(protPoolStoichiometry), ' ECM_protein_pool_biomass[e]', ...
%       ' => ECM_biomass[e]'];
   
rxnsToAdd.rxns = {'ECM_biomass'};
rxnsToAdd.equations = {rxn};
model2 = addRxns(model2,rxnsToAdd, 3);

%And finally, the total biomass function
%So, we still strive for that the biomass equation should have the unit
%growth per hour, so the stoichiometric constant for tumor growth should be
%1, while the ECM should be adapted. That should be
%fracECM/((1-fracECM)*fracC)
sc = fracECM/((1-fracECM)*fracC);

rxn = ['biomass[c] + ', ...
       num2str(sc), ' ECM_biomass[f_e]', ...
       ' =>'];
rxnsToAdd.rxns = {'total_biomass'};
rxnsToAdd.equations = {rxn};
model2 = addRxns(model2,rxnsToAdd, 3);


outModel = model2;

%for test:
%ecModelCorr = ltModel;
%fracC = 0.6;
%fracF = 0.2;
%fracECM = 0.2;
%fracGAGsInECM = 0.2;
