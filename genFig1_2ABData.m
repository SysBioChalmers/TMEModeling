cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling

%load the "small" model (i.e. just the ltModel, one cell type only)
load('data/ltModel.mat');

cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"

bloodData = prepBloodData();

%get the metabolites for each exchange reaction
[~, exchRxnIndAll] = getExchangeRxns(ltModel, 'in');
exchRxnMetsAll = cell(length(exchRxnIndAll), 1);

for i = 1:length(exchRxnMetsAll)
   exchRxnMetsAll{i} = ltModel.metNames{ltModel.S(:, exchRxnIndAll(i)) == 1};
end

%filter mets and exchange rxns that does not exist in the blood data (it is meaningless to work with them)
%filters 5 metabolites + prot_pool
sel = ismember(strcat(exchRxnMetsAll,'[e]'), bloodData.totMets);
exchRxnMetsAll(~sel)%H2O, Pi, sulfate, Fe2+, hypoxanthine, prot_pool
%check the other way around - all exist there
%sel2 = ismember(bloodData.totMets, strcat(exchRxnMets,'[s]'));
%bloodData.totMets(~sel2)%this is just water etc.



exchRxnMets = exchRxnMetsAll(sel);
exchRxnInd = exchRxnIndAll(sel);

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


%The "standard" range
a = (0.000001:0.000001:0.0001);
tic
D1_1 = runASimulation(ltModel, a, bloodData, cell_maintenance);
toc
%With 10^-9: Elapsed time is 174.101684 seconds.
%With 10^-6: Elapsed time is 199.900651 seconds.
save('data/D1_1.mat', 'D1_1')

aFigA = (0.000001:0.0000015:0.00015);
D2_0 = runASimulation(ltModel, aFigA, bloodData, cell_maintenance);
save('data/D2_0.mat', 'D2_0')

%rxnsToAdd = struct();
%rxnsToAdd.rxns = {'dho_exp'};
%rxnsToAdd.equations = {'(S)-dihydroorotate[c] =>'};
%ltModelWithDho = addRxns(ltModel,rxnsToAdd, 3);

%D1_1_DHO = runASimulation(ltModelWithDho, a, bloodData, cell_maintenance);
%save('data/D1_1_DHO.mat', 'D1_1_DHO')




D1_2 = runMetaboliteImportanceSimulation(ltModel, a, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, cell_maintenance, true);
save('data/D1_2.mat', 'D1_2');


%run metabolite importance simulation with a model with only 50% ATP costs:
%modify the biomass, NGAM and proteome constraint
ltModelMod = ltModel;
ngamVal = 1.833 * 0.5;
ltModelMod.lb(strcmp(ltModelMod.rxns,'MAR03964')) = ngamVal; %non-growth associated cell maintenance
ltModelMod.ub(strcmp(ltModelMod.rxns,'prot_pool_exchange')) = 8.775209e-03;%enzyme usage, a value fitted to the new biomass reaction and maintenance
%Do NOT remove and then add reactions, it will change the order! instead, alter the S matrix
ltModelMod.S(ltModelMod.S(:,strcmp(ltModelMod.rxns, 'MAR13082')) == 45, strcmp(ltModelMod.rxns, 'MAR13082')) = 45*0.5;
ltModelMod.S(ltModelMod.S(:,strcmp(ltModelMod.rxns, 'MAR13082')) == -45, strcmp(ltModelMod.rxns, 'MAR13082')) = -45*0.5;
constructEquations(ltModelMod, 'MAR13082')%ok
%check that we only changed this
%sum(sum(ltModelMod.S ~= ltModel.S,1),2)%5 items in the S matrix has changed, ok
D1_2C = runASimulation(ltModelMod, a, bloodData, ngamVal);
save('data/D1_2C.mat', 'D1_2C')
D1_2B = runMetaboliteImportanceSimulation(ltModelMod, a, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, ngamVal, true);
save('data/D1_2B.mat', 'D1_2B');

%Again, but 25%:
%modify the biomass, NGAM and proteome constraint
ltModelMod = ltModel;
ngamVal = 1.833 * 0.25;
ltModelMod.lb(strcmp(ltModelMod.rxns,'MAR03964')) = ngamVal; %non-growth associated cell maintenance
ltModelMod.ub(strcmp(ltModelMod.rxns,'prot_pool_exchange')) = 4.611405e-03;%enzyme usage, a value fitted to the new biomass reaction and maintenance
%Do NOT remove and then add reactions, it will change the order! instead, alter the S matrix
ltModelMod.S(ltModelMod.S(:,strcmp(ltModelMod.rxns, 'MAR13082')) == 45, strcmp(ltModelMod.rxns, 'MAR13082')) = 45*0.25;
ltModelMod.S(ltModelMod.S(:,strcmp(ltModelMod.rxns, 'MAR13082')) == -45, strcmp(ltModelMod.rxns, 'MAR13082')) = -45*0.25;
constructEquations(ltModelMod, 'MAR13082')%ok
%check that we only changed this
%sum(sum(ltModelMod.S ~= ltModel.S,1),2)%5 items in the S matrix has changed, ok
D1_2E = runASimulation(ltModelMod, a, bloodData, ngamVal);
save('data/D1_2E.mat', 'D1_2E')
D1_2D = runMetaboliteImportanceSimulation(ltModelMod, a, bloodData, exchRxnMets, mappingExchMets, exchRxnInd, ngamVal, true);
save('data/D1_2D.mat', 'D1_2D');






D1_3 = runASimulation(ltModel, a, bloodData, cell_maintenance, true);
save('data/D1_3.mat', 'D1_3')


%experiment with limiting lactate output, to keep pH in check
glucoseInd = find(strcmp(bloodData.totMets, 'glucose[e]'));
D1_8 = runASimulation(ltModel, a, bloodData, cell_maintenance, false, bloodData.totDxC(glucoseInd)/2);
save('data/D1_8.mat', 'D1_8')


%experiment with removing things from the biomass function to get an idea about what is most limiting - look at how the growth rate changes
%full biomass function:
%'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'

%remove ATP
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_9 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_9.mat', 'D1_9');

%remove nucleotides
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 45 H2O[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_10 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_10.mat', 'D1_10');

%remove glycogen
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_11 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_11.mat', 'D1_11');

%remove cofactors
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_12 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_12.mat', 'D1_12');

%remove protein
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_13 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_13.mat', 'D1_13');

%remove lipids
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_14 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_14.mat', 'D1_14');

%remove metabolites
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_15 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_15.mat', 'D1_15');

%now a special case, where we remove only the ATP cost of the protein
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass2[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_16 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_16.mat', 'D1_16');

%now a special case 2, where we remove the ATP cost of the protein + the other ATP cost
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass2[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_17 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_17.mat', 'D1_17');

%now a special case 3, where we remove the ATP cost of the protein + the other ATP cost + protein
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_18 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_18.mat', 'D1_18');


%now a special case 7, where we remove the ATP cost of the protein + the other ATP cost + lipids are removed
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass2[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_22 = runASimulation(tmpModel, a, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_22.mat', 'D1_22');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations using the blood flow model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling

%load the "small" model (i.e. just the ltModel, one cell type only)
load('data/ltModel.mat');

cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"

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


%The "standard" range is roughly 10 times higher here
a_bf = (0.00001:0.00001:0.001);
tic
D1_1_bloodflow = runASimulation(ltModel, a_bf, bloodData2, cell_maintenance);
toc
save('data/D1_1_bloodflow.mat', 'D1_1_bloodflow')
%The code below is for understanding why the spike in lipid uptake occurs
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'acetyl-CoA', false, 10^-3);%All comes from lipids
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'acetyl-CoA', true, 10^-3);%acetyl-CoA goes to citrate, i.e., TCA cycle
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'citrate', true, 10^-3);%most goes to isocitrate[m], 7% exported to cytosol where it is turned to isocitrate[c]
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'isocitrate', true, 10^-3);%turned to AKG in [m] and [c], produces NADPH
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'AKG', false, 10^-3);%Most comes from lipids
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'AKG', true, 10^-3);%most is turned into succinyl-CoA, i.e., TCA cycle
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'succinyl-CoA', true, 10^-3);%TCA cycle, ATP generation
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'NADH', true, 10^-3);%a small flux through complex 1 and something with DHAP, the large fluxes go through PYCR and lactate dehydrogenase
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'NADPH', true, 10^-3);% 4 things - major ones are:
%* 0.17168 MAR00478 DHAP[c] + H+[c] + NADPH[c] + 4.3481e-05 prot_pool[c] => NADP+[c] + sn-glycerol-3-phosphate[c]
%* 0.16037 MAR04335 dihydrofolate[m] + H+[m] + NADPH[m] + 0.00011684 prot_pool[c] => NADP+[m] + THF[m]
%* 0.04135 MAR00781 H+[c] + hexadecenal[c] + NADPH[c] + 0.00027666 prot_pool[c] => hexadecanal[c] + NADP+[c]
%* 0.00313 5,10-methenyl-THF[c] + NADPH[c] + 0.00021053 prot_pool[c] => 5,10-methylene-THF[c] + NADP+[c]
%These may represent alternative pathways to get rid of electrons
%Investigate DHAP first
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'sn-glycerol-3-phosphate', true, 10^-3);% forms a loop to go from NADPH to ubiquinol, i.e., similar to complex I, but no H+[i] gain, but cheaper from a protein pool perspective
%Now dihydrofolate
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'dihydrofolate', false, 10^-3);% 
listMetRxnsWithFluxes(ltModel, D1_1_bloodflow.resultSolutions{5}, 'THF', false, 10^-3);% 
%There is a loop with dihydrofolate and THF that in practice transports NADPH to the cytosol, to be turned into ubiquinol by the DHAP loop, thereby bypassing complex I.
%Conclusion: Lipids provide a higher FADH2/NADH ratio, and are therefore more beneficial if we want to bypass complex I, since less effort needs to be spent 
%on getting rid of the NADH.

D1_2_bloodflow = runMetaboliteImportanceSimulation(ltModel, a_bf, bloodData2, exchRxnMets, mappingExchMets, exchRxnInd, cell_maintenance, true);
save('data/D1_2_bloodflow.mat', 'D1_2_bloodflow');


rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_9_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_9_bloodflow.mat', 'D1_9_bloodflow');

%remove nucleotides
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 45 H2O[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_10_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_10_bloodflow.mat', 'D1_10_bloodflow');

%remove glycogen
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_11_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_11_bloodflow.mat', 'D1_11_bloodflow');

%remove cofactors
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_12_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_12_bloodflow.mat', 'D1_12_bloodflow');

%remove protein
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_13_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_13_bloodflow.mat', 'D1_13_bloodflow');

%remove lipids
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_14_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_14_bloodflow.mat', 'D1_14_bloodflow');

%remove metabolites
rxnsToAdd.rxns = {'incomplete_biomass'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]'};
tmpModel = addRxns(ltModel,rxnsToAdd, 3);
D1_15_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_15_bloodflow.mat', 'D1_15_bloodflow');

%now a special case, where we remove only the ATP cost of the protein
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'45 ATP[c] + 0.0267 DNA[n] + 45 H2O[c] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass2[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => 45 ADP[c] + 45 H+[c] + 45 Pi[c] + biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_16_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_16_bloodflow.mat', 'D1_16_bloodflow');

%now a special case 2, where we remove the ATP cost of the protein + the other ATP cost
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass2[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_17_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_17_bloodflow.mat', 'D1_17_bloodflow');

%now a special case 3, where we remove the ATP cost of the protein + the other ATP cost + protein
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 0.2212 lipid_pool_biomass[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_18_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_18_bloodflow.mat', 'D1_18_bloodflow');


%now a special case 7, where we remove the ATP cost of the protein + the other ATP cost + lipids are removed
metsToAdd.mets = {'protein_pool_biomass2'};
metsToAdd.metNames = {'protein_pool_biomass2'};
metsToAdd.compartments = {'c'};
tmpModel = addMets(ltModel, metsToAdd, false);
rxnsToAdd.rxns = {'incomplete_biomass', 'special_protein_pool'};
rxnsToAdd.equations = {'0.0267 DNA[n] + 0.1124 RNA[c] + 0.4062 glycogen[c] + 0.0012 cofactor_pool_biomass[c] + 5.3375 protein_pool_biomass2[c] + 0.4835 metabolite_pool_biomass[c] => biomass[c]', ...
    '0.0721 glycine[c] + 0.0801 alanine[c] + 0.0512 arginine[c] + 0.0375 asparagine[c] + 0.0556 aspartate[c] + 0.0183 cysteine[c] + 0.0428 glutamine[c] + 0.0783 glutamate[c] + 0.0228 histidine[c] + 0.0442 isoleucine[c] + 0.0911 leucine[c] + 0.0719 lysine[c] + 0.0222 methionine[c] + 0.0368 phenylalanine[c] + 0.051 proline[c] + 0.0661 serine[c] + 0.0535 threonine[c] + 0.0098 tryptophan[c] + 0.0281 tyrosine[c] + 0.0667 valine[c] => protein_pool_biomass2[c]'};
tmpModel = addRxns(tmpModel,rxnsToAdd, 3);
D1_22_bloodflow = runASimulation(tmpModel, a_bf, bloodData, cell_maintenance, false, NaN, 'incomplete_biomass');
save('data/D1_22_bloodflow.mat', 'D1_22_bloodflow');




