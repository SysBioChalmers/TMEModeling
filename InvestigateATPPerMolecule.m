%Generates data for table 1 and associated supplementary figures.
cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling

%load the models
load('data/ltModel.mat');

substrateNames = { ...
                  'lactate', ...
                  'aspartate','glutamine','glycine','proline','serine','threonine',...
                  'alanine','arginine','asparagine','cysteine','glutamate','histidine', ...
                  'isoleucine','leucine','lysine','methionine','phenylalanine','tryptophan','tyrosine','valine'
};
substrateRxns = { ...
                  'MAR09135_REV', ... %lactate
                  'MAR09070_REV', ... %aspartate
                  'MAR09063_REV', ... %glutamine
                  'MAR09067_REV', ... %glycine
                  'MAR09068_REV', ... %proline
                  'MAR09069_REV', ... %serine
                  'MAR09044_REV', ... %threonine
                  'MAR09061_REV', ... %alanine
                  'MAR09066_REV', ... %arginine
                  'MAR09062_REV', ... %asparagine
                  'MAR09065_REV', ... %cysteine
                  'MAR09071_REV', ... %glutamate
                  'MAR09038_REV', ... %histidine
                  'MAR09039_REV', ... %isoleucine
                  'MAR09040_REV', ... %leucine
                  'MAR09041_REV', ... %lysine
                  'MAR09042_REV', ... %methionine
                  'MAR09043_REV', ... %phenylalanine
                  'MAR09045_REV', ... %tryptophan
                  'MAR09064_REV', ... %tyrosine
                  'MAR09046_REV'  ... %valine
};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First look at ATP production per substrate in hypoxia without enzyme constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mBase = ltModel;


%optimize for ATP, test to send in different substrates
%first close all metabolite uptake reactions
exchRxns = getExchangeRxns(mBase, 'in');
exchRxnsInd = ismember(mBase.rxns, exchRxns);
mBase.ub(exchRxnsInd) = 0;

%now, optimize for ATP, when supplying lactate and glutamine, with limited oxygen and see what it prefers!
%I reason lactate is similar to glucose, except that we don't get any extra ATP from glycolysis, 
%but redox should be similar.
mBase.lb(strcmp(ltModel.rxns,'MAR03964')) = 0;%remove NGAM
mBase.ub(length(mBase.ub)) = Inf; %disable enzyme constraints
mBase.ub(strcmp(mBase.rxns, 'MAR09047_REV')) = Inf;%H2O
mBase.ub(strcmp(mBase.rxns, 'MAR09072_REV')) = Inf;%Pi
mBase.ub(strcmp(mBase.rxns, 'MAR09079_REV')) = Inf;%H+
mBase.ub(strcmp(mBase.rxns, 'MAR09048_REV')) = 1;%oxygen
mBase.ub(strcmp(mBase.rxns, 'MAR09034_REV')) = 0;%glucose
mBase.c = double(strcmp(mBase.rxns,'MAR03964'));%set objective function to ATP consumption

nsub = length(substrateRxns);
resHypoxiaWithProdh = NaN(nsub,1);
resHypoxiaNoProdh = NaN(nsub,1);
resHypoxiaNoProdhNoO2 = NaN(nsub,1);
sols = cell(nsub,1);

for i = 1:nsub
    disp(i)
    m = mBase;
    m.ub(strcmp(mBase.rxns, substrateRxns{i})) = 10;%the substrate
    res = solveLP(m,1);
    resHypoxiaWithProdh(i) = -res.f;%ATP prod
    %now block the reversed prodh reaction
    m.ub(strcmp(mBase.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model
    res = solveLP(m,1);
    resHypoxiaNoProdh(i) = -res.f;%ATP prod
    sols{i} = res;
    m.ub(strcmp(mBase.rxns, 'MAR09048_REV')) = 0;%Now also block oxygen
    res = solveLP(m,1);
    resHypoxiaNoProdhNoO2(i) = -res.f;%ATP prod
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then look at ATP production per substrate with enzyme constraints, but no hypoxia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mBase = ltModel;

%optimize for ATP, test to send in different substrates
%first close all metabolite uptake reactions
exchRxns = getExchangeRxns(mBase, 'in');
exchRxnsInd = ismember(mBase.rxns, exchRxns);
mBase.ub(exchRxnsInd) = 0;

%now, optimize for ATP, when supplying lactate and glutamine, with limited oxygen and see what it prefers!
%I reason lactate is similar to glucose, except that we don't get any extra ATP from glycolysis, 
%but redox should be similar.

mBase.lb(strcmp(ltModel.rxns,'MAR03964')) = 0;%remove NGAM
mBase.ub(length(mBase.ub)) = 0.001; %enable enzyme constraints at an appropriate level
mBase.ub(strcmp(mBase.rxns, 'MAR09047_REV')) = Inf;%H2O
mBase.ub(strcmp(mBase.rxns, 'MAR09072_REV')) = Inf;%Pi
mBase.ub(strcmp(mBase.rxns, 'MAR09079_REV')) = Inf;%H+
mBase.ub(strcmp(mBase.rxns, 'MAR09048_REV')) = Inf;%oxygen
mBase.ub(strcmp(mBase.rxns, 'MAR09034_REV')) = 0;%glucose
mBase.c = double(strcmp(mBase.rxns,'MAR03964'));%set objective function to ATP consumption

nsub = length(substrateRxns);
resEnzWithProdh = NaN(nsub,1);
resEnzNoProdh = NaN(nsub,1);
sols2 = cell(nsub,1);

for i = 1:nsub
    disp(i)
    m = mBase;
    m.ub(strcmp(mBase.rxns, substrateRxns{i})) = 10;%the substrate
    res = solveLP(m,1);
    sols2{i} = res;
    resEnzWithProdh(i) = -res.f;%ATP prod
    %now block the reversed prodh reaction
    m.ub(strcmp(mBase.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model
    res = solveLP(m,1);
    resEnzNoProdh(i) = -res.f;%ATP prod
end

table(substrateNames.', resHypoxiaWithProdh, resHypoxiaNoProdh, resHypoxiaNoProdhNoO2, resEnzWithProdh, resEnzNoProdh)


%create a list with TCA cycle fluxes and percentage of flux through complex I vs complex III
s = struct();
s.TCACycleFluxesLowO2NoProdh = NaN(nsub,1);
s.ATPLowO2NoProdh = resHypoxiaNoProdh;
s.complexILowO2NoProdh = NaN(nsub,1);
s.complexIIILowO2NoProdh = NaN(nsub,1);
s.complexVLowO2NoProdh = NaN(nsub,1);
s.TCACycleFluxesEnzLim = NaN(nsub,1);
s.complexIEnzLim = NaN(nsub,1);
s.complexIIIEnzLim = NaN(nsub,1);
s.complexVEnzLim = NaN(nsub,1);
s.ATPEnzLim = resEnzNoProdh;
s.mets = substrateNames;

for i = 1:nsub
    TCAGTP = sols{i}.x(strcmp(m.rxns, 'MAR04147_REV'));
    TCAATP = sols{i}.x(strcmp(m.rxns, 'MAR04152'));
    s.TCACycleFluxesLowO2NoProdh(i) = TCAGTP + TCAATP;
    s.complexILowO2NoProdh(i) = sols{i}.x(strcmp(m.rxns, 'MAR06921'));
    s.complexIIILowO2NoProdh(i) = sols{i}.x(strcmp(m.rxns, 'MAR06918'));
    s.complexVLowO2NoProdh(i) = sols{i}.x(strcmp(m.rxns, 'MAR06916'));

    TCAGTP = sols2{i}.x(strcmp(m.rxns, 'MAR04147_REV'));
    TCAATP = sols2{i}.x(strcmp(m.rxns, 'MAR04152'));
    s.TCACycleFluxesEnzLim(i) = TCAGTP + TCAATP;
    s.complexIEnzLim(i) = sols2{i}.x(strcmp(m.rxns, 'MAR06921'));
    s.complexIIIEnzLim(i) = sols2{i}.x(strcmp(m.rxns, 'MAR06918'));
    s.complexVEnzLim(i) = sols2{i}.x(strcmp(m.rxns, 'MAR06916'));
end

save('data/D3_1.mat', 's')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce similar data for running complex II in reverse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 cases:
%1. Hypoxia, citrate synthase in reverse closed
%2. Hypoxia, citrate synthase in reverse open
%3. Enzyme constraints, citrate synthase in reverse closed
%4. Enzyme constraints, citrate synthase in reverse open
%We shut down PRODH in reverse in this case

%add succinate export
rxnsToAdd = struct();
rxnsToAdd.rxns = {'succExp'};
rxnsToAdd.equations = {'succinate[c] => succinate[e]'};
rxnsToAdd.ub = Inf;
mBase = addRxns(ltModel, rxnsToAdd, 3);
%shut down PRODH in reverse
mBase.ub(strcmp(mBase.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model

%optimize for ATP, test to send in different substrates
%first close all metabolite uptake reactions
exchRxns = getExchangeRxns(mBase, 'in');
exchRxnsInd = ismember(mBase.rxns, exchRxns);
mBase.ub(exchRxnsInd) = 0;

mBase.lb(strcmp(ltModel.rxns,'MAR03964')) = 0;%remove NGAM
mBase.ub(strcmp(mBase.rxns, 'MAR09047_REV')) = Inf;%H2O
mBase.ub(strcmp(mBase.rxns, 'MAR09072_REV')) = Inf;%Pi
mBase.ub(strcmp(mBase.rxns, 'MAR09079_REV')) = Inf;%H+
mBase.ub(strcmp(mBase.rxns, 'MAR09048_REV')) = 1;%oxygen
mBase.ub(strcmp(mBase.rxns, 'MAR09034_REV')) = 0;%glucose

mBase.c = double(strcmp(mBase.rxns,'MAR03964'));%set objective function to ATP consumption

baseModelSuccExp = mBase;

%First case 1 and 2:
%%%%%%%%%%%%%%%%%%%%

mBase2 = mBase;

mBase2.ub(strcmp(mBase.rxns, 'prot_pool_exchange')) = Inf; %disable enzyme constraints


nsub = length(substrateRxns);
resCIIHypoxiaNoCit = NaN(nsub,1);
resCIIHypoxia = NaN(nsub,1);
sols_c1 = cell(nsub,1);
sols_c2 = cell(nsub,1);

for i = 1:nsub
    disp(i)
    %Case 1 - block citrate synthase in reverse
    m = mBase2;
    m.ub(strcmp(mBase.rxns, substrateRxns{i})) = 10;%the substrate
    m.ub(strcmp(mBase.rxns, 'MAR04145')) = 0;%citrate synthase in reverse - the reaction is defined backwards
    sols_c1{i} = solveLP(m,1);
    resCIIHypoxiaNoCit(i) = -sols_c1{i}.f;%ATP prod
    %Case 2 - open citrate synthase in reverse
    m.ub(strcmp(mBase.rxns, 'MAR04145')) = Inf;%citrate synthase in reverse - the reaction is defined backwards
    sols_c2{i} = solveLP(m,1);
    resCIIHypoxia(i) = -sols_c2{i}.f;%ATP prod
end

%Now case 3 and 4, i.e., enzyme constraints:
%%%%%%%%%%%%%%%%%%%%

%now, optimize for ATP, when supplying lactate and glutamine, with limited oxygen and see what it prefers!
%I reason lactate is similar to glucose, except that we don't get any extra ATP from glycolysis, 
%but redox should be similar.
mBase2 = mBase;
mBase2.ub(strcmp(mBase.rxns, 'prot_pool_exchange')) = 0.001; %enable enzyme constraints at an appropriate level
mBase2.ub(strcmp(mBase.rxns, 'MAR09048_REV')) = Inf;%oxygen - no hypoxia here

nsub = length(substrateRxns);
resCIIEnzNoCit = NaN(nsub,1);
resCIIEnz = NaN(nsub,1);
sols_c3 = cell(nsub,1);
sols_c4 = cell(nsub,1);

for i = 1:nsub
    disp(i)
    m = mBase2;
    %Case 3 - block citrate synthase in reverse
    m.ub(strcmp(mBase.rxns, substrateRxns{i})) = 10;%the substrate
    m.ub(strcmp(mBase.rxns, 'MAR04145')) = 0;%citrate synthase in reverse - the reaction is defined backwards
    sols_c3{i} = solveLP(m,1);
    resCIIEnzNoCit(i) = -sols_c3{i}.f;%ATP prod
    %Case 4 - open citrate synthase in reverse
    m.ub(strcmp(mBase.rxns, 'MAR04145')) = Inf;%citrate synthase in reverse - the reaction is defined backwards
    sols_c4{i} = solveLP(m,1);
    resCIIEnz(i) = -sols_c4{i}.f;%ATP prod
end

table(substrateNames.', resCIIHypoxiaNoCit, resCIIHypoxia, resCIIEnzNoCit, resCIIEnz)
resCIIHypoxiaNoCit - resCIIHypoxia %they are the same
resCIIEnzNoCit - resCIIEnz %they are the same

%Check what pathways are chosen for lactate
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{1}, 'L-lactate', true, 10^-4);
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{1}, 'pyruvate', true, 10^-4); %60% goes to OAA at the cost of ATP
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{1}, 'OAA', true, 10^-4); %to malate
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{1}, 'malate', true, 10^-4); %to fumarate
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{1}, 'fumarate', true, 10^-4); %to succinate, i.e., complex II in reverse 6
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{1}, 'succinate', true, 10^-4); %succinate exported
%So, we lose 1 ATP and gain 1.33 ATP since complex I pumps 4 protons, given that we have enough NADH.

%The same for glutamine
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{3}, 'glutamine', true, 10^-4); %imported and converted to glutamate
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{3}, 'glutamate', true, 10^-4); %10% -> AKG, 45% -> L-glutamate 5-semialdehyde[m], 45% används med L-glutamate 5-semialdehyde[m] för att bilda AKG + ornithine
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{3}, 'AKG', true, 10^-4); %goes to succinyl-CoA
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{3}, 'succinyl-CoA', true, 10^-4); %goes to succinate
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{3}, 'succinate', true, 10^-4); %succinate is exported - this is 55% of the glutamine
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{3}, 'ornithine', true, 10^-4); %is eventually exported
listMetRxnsWithFluxes(baseModelSuccExp, sols_c1{3}, 'citrulline', true, 10^-4); %is eventually exported


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now, investigate the PRODH and complex II in reverse cases with the diffusion model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We look at these cases:
%1. None of them on
%2. PRODH in reverse on
%3. Succinate export on, also resulting in that complex II can be run in reverse. Citrate synthase in reverse is blocked.
%4. Succinate export on, also resulting in that complex II can be run in reverse. Citrate synthase in reverse is open.
%5. All advantageous reactions are open.

%add succinate export
rxnsToAdd = struct();
rxnsToAdd.rxns = {'succExp'};
rxnsToAdd.equations = {'succinate[c] => succinate[e]'};
rxnsToAdd.ub = Inf;
mBase = addRxns(ltModel, rxnsToAdd, 3);


%Setup the diffusion model
%%%%%%%%%%%%%%%%%%%%%%%%%%

bloodData = prepBloodData();

%get the metabolites for each exchange reaction
[~, exchRxnIndAll] = getExchangeRxns(mBase, 'in');
exchRxnMetsAll = cell(length(exchRxnIndAll), 1);

for i = 1:length(exchRxnMetsAll)
   exchRxnMetsAll{i} = m.metNames{mBase.S(:, exchRxnIndAll(i)) == 1};
end

%filter mets and exchange rxns that does not exist in the blood data (it is meaningless to work with them)
%filters 5 metabolites + prot_pool
sel = ismember(strcat(exchRxnMetsAll,'[e]'), bloodData.totMets);
exchRxnMetsAll(~sel)%H2O, Pi, sulfate, Fe2+, hypoxanthine, prot_pool

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

%Set NGAM to zero
mBase.lb(strcmp(mBase.rxns,'MAR03964')) = 0;%no NGAM

%%%%%%%%%%%%%%%%%%
% case 1 - nothing on
%%%%%%%%%%%%%%%%%%

m = mBase;
%shut down PRODH in reverse
m.ub(strcmp(mBase.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model
%shut down succinate export (which stops running complex I in reverse)
m.ub(strcmp(mBase.rxns, 'succExp')) = 0;

tic
D5_6_1 = runASimulation(m, a, bloodData, 0, false, NaN, 'MAR03964'); 
toc
save('data/D5_6_1.mat', 'D5_6_1')

%%%%%%%%%%%%%%%%%%
% case 2 - PRODH in reverse on
%%%%%%%%%%%%%%%%%%

m = mBase;
%shut down succinate export (which stops running complex I in reverse)
m.ub(strcmp(mBase.rxns, 'succExp')) = 0;

tic
D5_6_2 = runASimulation(m, a, bloodData, 0, false, NaN, 'MAR03964'); 
toc
save('data/D5_6_2.mat', 'D5_6_2')

%%%%%%%%%%%%%%%%%%
% case 3 - Succinate export on, also resulting in that complex II can be run in reverse. Citrate synthase in reverse is blocked.
%          This is not used in the final figures.
%%%%%%%%%%%%%%%%%%

m = mBase;
%shut down PRODH in reverse
m.ub(strcmp(mBase.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model
%shut down citrate synthase in reverse
m.ub(strcmp(mBase.rxns, 'MAR04145')) = 0;%citrate synthase in reverse - the reaction is defined backwards

tic
D5_6_3 = runASimulation(m, a, bloodData, 0, false, NaN, 'MAR03964'); 
toc
save('data/D5_6_3.mat', 'D5_6_3')

%%%%%%%%%%%%%%%%%%
% case 4 - Succinate export on, also resulting in that complex II can be run in reverse. Citrate synthase in reverse is open.
%%%%%%%%%%%%%%%%%%

m = mBase;
%shut down PRODH in reverse
m.ub(strcmp(mBase.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model

tic
D5_6_4 = runASimulation(m, a, bloodData, 0, false, NaN, 'MAR03964'); 
toc
save('data/D5_6_4.mat', 'D5_6_4')

listMetRxnsWithFluxes(baseModelSuccExp, D5_6_4.resultSolutions{60}, 'succinate', false, 10^-4); %succinate is exported - this is 55% of the glutamine


%%%%%%%%%%%%%%%%%%
% case 5 - All open
%%%%%%%%%%%%%%%%%%

m = mBase;

tic
D5_6_5 = runASimulation(m, a, bloodData, 0, false, NaN, 'MAR03964'); 
toc
save('data/D5_6_5.mat', 'D5_6_5')


