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
% First look at glutamine addiction in hypoxia without enzyme constraints
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

mBase2 = convertToIrrev(ihuman);

%test galactose, temp
exchRxns2 = getExchangeRxns(mBase2, 'in');
exchRxnsInd2 = ismember(mBase2.rxns, exchRxns2);
mBase2.ub(exchRxnsInd2) = 0;

%now, optimize for ATP, when supplying lactate and glutamine, with limited oxygen and see what it prefers!
%I reason lactate is similar to glucose, except that we don't get any extra ATP from glycolysis, 
%but redox should be similar.
%mBase2.lb(strcmp(ltModel.rxns,'MAR03964')) = 0;%remove NGAM
%mBase2.ub(length(mBase.ub)) = Inf; %disable enzyme constraints
mBase2.ub(strcmp(mBase2.rxns, 'MAR09047_REV')) = Inf;%H2O
mBase2.ub(strcmp(mBase2.rxns, 'MAR09072_REV')) = Inf;%Pi
mBase2.ub(strcmp(mBase2.rxns, 'MAR09079_REV')) = Inf;%H+
mBase2.ub(strcmp(mBase2.rxns, 'MAR09048_REV')) = 1;%oxygen
mBase2.ub(strcmp(mBase2.rxns, 'MAR09034_REV')) = 0;%glucose
mBase2.c = double(strcmp(mBase2.rxns,'MAR03964'));%set objective function to ATP consumption


m = mBase2;
m.ub(strcmp(mBase2.rxns, 'MAR09140_REV')) = 10;%galactose
%m.ub(strcmp(mBase2.rxns, 'MAR09034_REV')) = 10;%glucose
m.ub(strcmp(mBase2.rxns, 'MAR09048_REV')) = 0;%oxygen
%m.c = double(strcmp(mBase2.rxns,'MAR04130'));%set objective function to ATP consumption
%m.c = double(strcmp(mBase.rxns,'MAR09140_REV'));%set objective function to ATP consumption

%constructEquations(m, 'MAR09140')

res = solveLP(m,1);
-res.f


table(substrateNames.', resHypoxiaWithProdh, resHypoxiaNoProdh,resHypoxiaNoProdhNoO2)

%investigate why glycine is beneficial
listMetRxnsWithFluxes(mBase, sols{1}, 'acetyl-CoA', true, 10^-9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then look at glutamine addiction with enzyme constraints, but no hypoxia
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
%Now investigate the role of dihydroorotate. This looks very interesting. It may be the 
%same thing as PRODH, that it is possible to run complex I without oxygen.

mBase2 = convertToIrrev(ihuman);


%optimize for ATP, test to send in different substrates
%first close all metabolite uptake reactions
exchRxns = getExchangeRxns(mBase2, 'in');
exchRxnsInd = ismember(mBase2.rxns, exchRxns);
mBase2.ub(exchRxnsInd) = 0;

%now, optimize for ATP, when supplying lactate and glutamine, with limited oxygen and see what it prefers!
%I reason lactate is similar to glucose, except that we don't get any extra ATP from glycolysis, 
%but redox should be similar.

%mBase2.lb(strcmp(ltModel.rxns,'MAR03964')) = 0;%remove NGAM
%mBase2.ub(length(mBase2.ub)) = Inf; %disable enzyme constraints
mBase2.ub(strcmp(mBase2.rxns, 'MAR09047_REV')) = Inf;%H2O
mBase2.ub(strcmp(mBase2.rxns, 'MAR09072_REV')) = Inf;%Pi
mBase2.ub(strcmp(mBase2.rxns, 'MAR09079_REV')) = Inf;%H+
mBase2.ub(strcmp(mBase2.rxns, 'MAR09048_REV')) = 1;%oxygen
mBase2.ub(strcmp(mBase2.rxns, 'MAR09034_REV')) = 0;%glucose
mBase2.c = double(strcmp(mBase2.rxns,'MAR03964'));%set objective function to ATP consumption

%add export of (S)-dihydroorotate
rxnsToAdd = struct();
rxnsToAdd.rxns = {'dho_exp'};
rxnsToAdd.equations = {'(S)-dihydroorotate[c] =>'};
mBase3 = addRxns(mBase2,rxnsToAdd, 3);


nsub = length(substrateRxns);
resHypoxiaWithProdh3 = NaN(nsub,1);
resHypoxiaNoProdh3 = NaN(nsub,1);
resHypoxiaNoProdhNoO23 = NaN(nsub,1);
resHypoxiaNoProdhNoO2NoDHO3 = NaN(nsub,1);
sols3 = cell(nsub,1);

for i = 1:nsub
    disp(i)
    m = mBase2;
    m.ub(strcmp(m.rxns, substrateRxns{i})) = 10;%the substrate
    res = solveLP(m,1);
    resHypoxiaWithProdh3(i) = -res.f;%ATP prod
    %now block the reversed prodh reaction
    m.ub(strcmp(m.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model
    res = solveLP(m,1);
    resHypoxiaNoProdh3(i) = -res.f;%ATP prod
    m.ub(strcmp(mBase.rxns, 'MAR09048_REV')) = 0;%Now also block oxygen
    res = solveLP(m,1);
    resHypoxiaNoProdhNoO23(i) = -res.f;%ATP prod
    
    m = mBase3;
    m.ub(strcmp(m.rxns, substrateRxns{i})) = 10;%the substrate
    %now block the reversed prodh reaction
    m.ub(strcmp(m.rxns, 'MAR03838')) = 0;%so, the prodh reaction is actually defined in reverse in the model
    res = solveLP(m,1);
    resHypoxiaNoProdhNoO2NoDHO3(i) = -res.f;%ATP prod
    sols3{i} = res;
end

table(substrateNames.', resHypoxiaWithProdh3, resHypoxiaNoProdh3,resHypoxiaNoProdhNoO23, resHypoxiaNoProdhNoO2NoDHO3)

listMetRxnsWithFluxes(mBase3, sols3{3}, '(S)-dihydroorotate', true, 10^-9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%below is code for understanding what is happening in the models


listMetRxnsWithFluxes(mBase, sols{1}, 'L-lactate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'pyruvate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'acetyl-CoA', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'acetate', true, 10^-9);

listMetRxnsWithFluxes(mBase, sols{3}, 'glutamate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{3}, 'AKG', true, 10^-9);



listMetRxnsWithFluxes(mBase, sols2{1}, 'L-lactate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols2{1}, 'pyruvate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols2{1}, 'acetyl-CoA', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols2{1}, 'citrate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols2{1}, 'isocitrate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols2{1}, 'AKG', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols2{1}, 'succinyl-CoA', true, 10^-9);%so, all flux goes through here
listMetRxnsWithFluxes(mBase, sols2{1}, 'NADH', true, 10^-9);%interestingly, pycr is also used for lactate!
listMetRxnsWithFluxes(mBase, sols2{1}, '1-pyrroline-5-carboxylate', false, 10^-9);%And, it forms a loop with PRODH to turn NADH into ubiquinol, thereby bypassing complex I!


listMetRxnsWithFluxes(mBase, sols2{3}, 'NADH', true, 10^-9);


%2 interesting things:
%1. Serine and threonine give high ATP yield in hypoxia (3.3333 extra)
%2: There is an odd behavior with serine, threonine, and glycine - they gain 0.66667 ATP extra in hypoxia without prodh

%1
%serine first
listMetRxnsWithFluxes(mBase, sols{6}, 'serine', true, 10^-9);
sols{6}.x(strcmp(m.rxns,'MAR04152')) %0.0096
%    {'MAR04198'}           4           {'pyruvate[c] + serine[c] + 0.00041197 prot_pool[c] => alanine[c] + hydroxypyruvate[c]'}
%    {'MAR00622'}           6           {'PE-LD pool[c] + serine[c] + 0.010081 prot_pool[c] => ethanolamine[c] + PS-LD pool[c]'}
%    {'MAR05079'}           6           {'serine[s] => serine[c]'                                                              }
%    {'MAR05474'}           4           {'alanine[c] + serine[s] + 0.00027751 prot_pool[c] => alanine[s] + serine[c]'          }
%So, some alanine is generated and discarded - hydroxypyruvate is more interesting
listMetRxnsWithFluxes(mBase, sols{6}, 'hydroxypyruvate', true, 10^-9);
%{'MAR08782'}           4           {'H+[c] + hydroxypyruvate[c] + NADH[c] + 1.9427e-05 prot_pool[c] => glycerate[c] + NAD+[c]'}
%Gets rid of NADH
listMetRxnsWithFluxes(mBase, sols{6}, 'glycerate', true, 10^-9);
% {'MAR08781'}           4           {'ATP[c] + glycerate[c] + 0.00027751 prot_pool[c] => 2-phospho-D-glycerate[c] + ADP[c] + H+[c]'}
%hmm, costs ATP
listMetRxnsWithFluxes(mBase, sols{6}, '2-phospho-D-glycerate', true, 10^-9);
%{'MAR04363'}           8           {'2-phospho-D-glycerate[c] + 0.00015979 prot_pool[c] => H2O[c] + PEP[c]'}
%we get PEP, also from another path
listMetRxnsWithFluxes(mBase, sols{6}, 'PEP', true, 10^-9);
%{'MAR04358'}           4           {'ADP[c] + H+[c] + PEP[c] + 5.5866e-06 prot_pool[c] => ATP[c] + pyruvate[c]'}
%{'MAR04210'}           4           {'H+[c] + PEP[c] + UDP[c] + 1.3615e-05 prot_pool[c] => pyruvate[c] + UTP[c]'}
%restores lost ATP, the other path unclear
listMetRxnsWithFluxes(mBase, sols{6}, 'pyruvate', true, 10^-9);
%{'MAR04388'}           4           {'H+[c] + NADH[c] + pyruvate[c] + 3.6479e-05 prot_pool[c] => L-lactate[c] + NAD+[c]'   }
%{'MAR04198'}           4           {'pyruvate[c] + serine[c] + 0.00041197 prot_pool[c] => alanine[c] + hydroxypyruvate[c]'}
%So, half loops back to the first reaction, the rest is exported as lactate to remove NADH
%Now, the ethanolamine[c]
listMetRxnsWithFluxes(mBase, sols{6}, 'ethanolamine', true, 10^-9);
%{'MAR00647'}           6           {'ethanolamine[c] + 0.00027751 prot_pool[c] => acetaldehyde[c] + H+[c] + NH3[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'acetaldehyde', true, 10^-9);
%{'MAR04398_REV'}           6           {'acetaldehyde[c] + GAP[c] + 0.00016392 prot_pool[c] => 2-deoxy-D-ribose-5-phosphate[c]'}
%consumes GAP here, where does that come from?
listMetRxnsWithFluxes(mBase, sols{6}, 'GAP', false, 10^-9);
%{'MAR04391'    }           4           {'DHAP[c] + 9.865e-07 prot_pool[c] => GAP[c]' - 1, part of glycolysis
%{'MAR04404'    }           2           {'D-xylulose-5-phosphate[c] + erythrose-4-phosphate[c] + 1.8904e-06 prot_pool[c] => fructose-6-phosphate[c] + GAP[c]'  } - 2
%{'MAR04501'    }           2           {'D-xylulose-5-phosphate[c] + ribose-5-phosphate[c] + 1.8904e-06 prot_pool[c] => GAP[c] + sedoheptulose-7-phosphate[c]'} - 3
%{'MAR04375_REV'}           4           {'fructose-1,6-bisphosphate[c] + 0.00018342 prot_pool[c] => DHAP[c] + GAP[c]'     -4, part of glycolysis
listMetRxnsWithFluxes(mBase, sols{6}, 'fructose-1,6-bisphosphate', false, 10^-9);
%{'MAR04301'}           4           {'fructose-6-phosphate[c] + UTP[c] + 2.2051e-05 prot_pool[c] => fructose-1,6-bisphosphate[c] + H+[c] + UDP[c]'}
%UDP cost, cancels out the one produced earlier
listMetRxnsWithFluxes(mBase, sols{6}, 'fructose-6-phosphate', false, 10^-9);
%{'MAR04404'}           2           {'D-xylulose-5-phosphate[c] + erythrose-4-phosphate[c] + 1.8904e-06 prot_pool[c] => fructose-6-phosphate[c] + GAP[c]'   } - same as 2 above
%{'MAR04565'}           2           {'GAP[c] + sedoheptulose-7-phosphate[c] + 1.0862e-06 prot_pool[c] => erythrose-4-phosphate[c] + fructose-6-phosphate[c]'} - loops back
listMetRxnsWithFluxes(mBase, sols{6}, 'D-xylulose-5-phosphate', false, 10^-9);
%{'MAR08992'}           4           {'2 ribulose-5-phosphate[c] + 5.7703e-07 prot_pool[c] => 2 D-xylulose-5-phosphate[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'ribulose-5-phosphate', false, 10^-9);
%{'MAR04352'}           4           {'ribose-5-phosphate[c] + 2.7752e-06 prot_pool[c] => ribulose-5-phosphate[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'ribose-5-phosphate', false, 10^-9);
%{'MAR04354'}           6           {'ribose-1-phosphate[c] + 0.0038793 prot_pool[c] => ribose-5-phosphate[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'ribose-1-phosphate', false, 10^-9);
%{'MAR04651'}           6           {'guanosine[c] + Pi[c] + 0.00028779 prot_pool[c] => guanine[c] + ribose-1-phosphate[c]'}
%So, here we get involved with guanosine again
listMetRxnsWithFluxes(mBase, sols{6}, '2-deoxy-D-ribose-5-phosphate', true, 10^-9);
%{'MAR04710_REV'}           6           {'2-deoxy-D-ribose-5-phosphate[c] + 0.00056966 prot_pool[c] => 2-deoxy-D-ribose-1-phosphate[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, '2-deoxy-D-ribose-1-phosphate', true, 10^-9);
%{'MAR04603'}           6           {'2-deoxy-D-ribose-1-phosphate[c] + guanine[c] + 0.00018587 prot_pool[c] => deoxyguanosine[c] + Pi[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'deoxyguanosine', true, 10^-9);
%{'MAR04601'}           6           {'ATP[m] + deoxyguanosine[m] + 0.0089043 prot_pool[c] => ADP[m] + dGMP[m] + H+[m]'}
%{'MAR04997'}           6           {'deoxyguanosine[c] => deoxyguanosine[m]'  
%So, a loss of ATP here
listMetRxnsWithFluxes(mBase, sols{6}, 'dGMP', true, 10^-9);
%{'MAR08481_REV'}           6           {'dATP[m] + dGMP[m] + 0.00027751 prot_pool[c] => dADP[m] + dGDP[m]'}
%Need to follow both of these, dADP first
listMetRxnsWithFluxes(mBase, sols{6}, 'dADP', true, 10^-9);
%{'MAR06615'}           6           {'dADP[m] + dTTP[m] + 0.00027751 prot_pool[c] => dATP[m] + dTDP[m]'}
%dATP loops back to one step above, follow dTDP
listMetRxnsWithFluxes(mBase, sols{6}, 'dTDP', true, 10^-9);
%{'MAR07886_REV'}           6           {'ATP[m] + dTDP[m] + 0.00027751 prot_pool[c] => ADP[m] + dTTP[m]'}
%Costs ATP
%loops back to a one step above, which in turn loops one step upward. The net effect is that dGMP is turned into
%dGDP at an ATP cost
%now follow the dGDP, which has now costed two ATP
listMetRxnsWithFluxes(mBase, sols{6}, 'dGDP', true, 10^-9);
%{'MAR02357'}           6           {'dGDP[m] + H2O[m] + oxidized thioredoxin[m] + 0.052061 prot_pool[c] => GDP[m] + thioredoxin[m]'}
%turned into GDP 
listMetRxnsWithFluxes(mBase, sols{6}, 'thioredoxin', true, 10^-9);
%{'MAR02358'}           6           {'NADP+[m] + thioredoxin[m] + 0.00033704 prot_pool[c] => H+[m] + NADPH[m] + oxidized thioredoxin[m]'}
%So, this loops back one step for generating NADPH. A bit strange, investigate more.
listMetRxnsWithFluxes(mBase, sols{6}, 'NADPH', true, 10^-9);
%{'MAR04335'}           2           {'dihydrofolate[m] + H+[m] + NADPH[m] + 0.00011684 prot_pool[c] => NADP+[m] + THF[m]'           }
%{'MAR04444'}           2           {'5,10-methenyl-THF[m] + NADPH[m] + 7.8555e-05 prot_pool[c] => 5,10-methylene-THF[m] + NADP+[m]'}
%{'MAR04656'}           2           {'folate[m] + NADPH[m] + 0.00011684 prot_pool[c] => dihydrofolate[m] + NADP+[m]'                }
%THe NADPH is restored in different ways. Follow these to see what is happening
listMetRxnsWithFluxes(mBase, sols{6}, 'THF', true, 10^-9);
%{'MAR04332_REV'}           2           {'NAD+[c] + THF[c] + 0.00011684 prot_pool[c] => dihydrofolate[c] + H+[c] + NADH[c]'}
%{'MAR03952_REV'}           2           {'THF[m] => THF[c]' 
%Generates NADH, loops back to above together with dihydro
listMetRxnsWithFluxes(mBase, sols{6}, '5,10-methylene-THF', true, 10^-9);
%{'MAR04440_REV'}           2           {'5,10-methylene-THF[m] + NAD+[m] + 2.0324e-05 prot_pool[c] => 5,10-methenyl-THF[m] + NADH[m]'}
%Generates NADH, loops back to above
listMetRxnsWithFluxes(mBase, sols{6}, 'dihydrofolate', true, 10^-9);
%{'MAR04440_REV'}           2           {'5,10-methylene-THF[m] + NAD+[m] + 2.0324e-05 prot_pool[c] => 5,10-methenyl-THF[m] + NADH[m]'}
%Generates NADH, loops back to above
%{'MAR04335'    }           2           {'dihydrofolate[m] + H+[m] + NADPH[m] + 0.00011684 prot_pool[c] => NADP+[m] + THF[m]'} - same as above
%{'MAR04654_REV'}           2           {'dihydrofolate[c] + NAD+[c] + 0.00011684 prot_pool[c] => folate[c] + NADH[c]'       }

%So, all these reactions generated NADH from ATP?

%different strategy, check if glycolysis is running
listMetRxnsWithFluxes(mBase, sols{6}, '1,3-bisphospho-D-glycerate', true, 10^-9);
%{'MAR04368'}           4           {'1,3-bisphospho-D-glycerate[c] + ADP[c] + 2.9149e-06 prot_pool[c] => 3-phospho-D-glycerate[c] + ATP[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'PEP', true, 10^-9);
%{'MAR04358'}           4           {'ADP[c] + H+[c] + PEP[c] + 5.5866e-06 prot_pool[c] => ATP[c] + pyruvate[c]'}
%{'MAR04210'}           4           {'H+[c] + PEP[c] + UDP[c] + 1.3615e-05 prot_pool[c] => pyruvate[c] + UTP[c]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'L-lactate', false, 10^-9);
%{'MAR04388'}           4           {'H+[c] + NADH[c] + pyruvate[c] + 3.6479e-05 prot_pool[c] => L-lactate[c] + NAD+[c]' }
%{'MAR06048'}           4           {'butyrate[s] + L-lactate[c] + 0.00027751 prot_pool[c] => butyrate[c] + L-lactate[s]'}
listMetRxnsWithFluxes(mBase, sols{6}, 'L-lactate', true, 10^-9)
%{'MAR06048'}           4           {'butyrate[s] + L-lactate[c] + 0.00027751 prot_pool[c] => butyrate[c] + L-lactate[s]'}
%{'MAR09135'}           4           {'L-lactate[s] => ' 
%So, glycolysis is running, some lactate is exported

%The conclusion is that this is a pretty massive network, but that it uses glycolysis to generate extra ATP

%Check if threonine is similar:
listMetRxnsWithFluxes(mBase, sols{7}, 'threonine', true, 10^-9);
%{'MAR04284'}           6           {'threonine[c] + 0.00027751 prot_pool[c] => acetaldehyde[c] + glycine[c]'          }
%{'MAR04348'}           4           {'threonine[c] + 3.3512e-05 prot_pool[c] => 2-oxobutyrate[c] + H+[c] + NH3[c]'     }
%{'MAR05470'}           6           {'glycine[c] + threonine[s] + 0.00027751 prot_pool[c] => glycine[s] + threonine[c]'}
%{'MAR05484'}           4           {'alanine[c] + threonine[s] + 0.00027751 prot_pool[c] => alanine[s] + threonine[c]'}
listMetRxnsWithFluxes(mBase, sols{7}, 'glycine', true, 10^-9);%exported
listMetRxnsWithFluxes(mBase, sols{7}, 'acetaldehyde', true, 10^-9);%exported
%{'MAR04398_REV'}           6           {'acetaldehyde[c] + GAP[c] + 0.00016392 prot_pool[c] => 2-deoxy-D-ribose-5-phosphate[c]'}
listMetRxnsWithFluxes(mBase, sols{7}, '2-deoxy-D-ribose-5-phosphate', true, 10^-9);
%{'MAR04710_REV'}           6           {'2-deoxy-D-ribose-5-phosphate[c] + 0.00056966 prot_pool[c] => 2-deoxy-D-ribose-1-phosphate[c]'}
listMetRxnsWithFluxes(mBase, sols{7}, '2-deoxy-D-ribose-1-phosphate', true, 10^-9);
%different strategy, check if glycolysis is running
listMetRxnsWithFluxes(mBase, sols{7}, '1,3-bisphospho-D-glycerate', true, 10^-9);
%{'MAR04368'}           4           {'1,3-bisphospho-D-glycerate[c] + ADP[c] + 2.9149e-06 prot_pool[c] => 3-phospho-D-glycerate[c] + ATP[c]'}
listMetRxnsWithFluxes(mBase, sols{7}, 'PEP', true, 10^-9);%exported
%{'MAR04210'}           4           {'H+[c] + PEP[c] + UDP[c] + 1.3615e-05 prot_pool[c] => pyruvate[c] + UTP[c]'}
listMetRxnsWithFluxes(mBase, sols{7}, 'L-lactate', true, 10^-9);%exported
%{'MAR06048'}        4.5e-06        {'butyrate[s] + L-lactate[c] + 0.00027751 prot_pool[c] => butyrate[c] + L-lactate[s]'}
%{'MAR09135'}        4.5e-06        {'L-lactate[s] => ' 

%So, threonine is also using glycolysis, but doesn't export much lactate, which is different from serine

%now, look at glycine, the second effect
listMetRxnsWithFluxes(mBase, sols{4}, 'glycine', true, 10^-9);
%generates half the amount of serine + extra NADH, so not surprising it can generate some extra ATP.

%Why is there a difference between glutamine and glutamate for the enzyme constrained case?
listMetRxnsWithFluxes(mBase, sols2{3}, 'glutamine', true, 10^-9);
%{'MAR04300'    }       0.092585        {'fructose-6-phosphate[c] + glutamine[c] + 1.7289e-05 prot_pool[c] => glucosamine-6-phosphate[c] + glutamate[c]'}
listMetRxnsWithFluxes(mBase, sols2{3}, 'glucosamine-6-phosphate', true, 10^-9);
%{'MAR04299'}       0.092585        {'glucosamine-6-phosphate[c] + H2O[c] + 3.321e-05 prot_pool[c] => fructose-6-phosphate[c] + H+[c] + NH3[c]'}
listMetRxnsWithFluxes(mBase, sols2{3}, 'NH3', true, 10^-9);%just exported
%look at glutamate
listMetRxnsWithFluxes(mBase, sols2{12}, 'glutamate', true, 10^-9);%just exported
%It is just about the protein cost of the exchange reactions, not very reliable information, they are probably the same

sols{1}.x(strcmp(m.rxns,'MAR04152')) %0.0096
sols{2}.x(strcmp(m.rxns,'MAR04152')) %0.0463

sols{1}.x(strcmp(m.rxns,'MAR04152')) - sols{2}.x(strcmp(m.rxns,'MAR04152')) %nothing

listMetRxnsWithFluxes(mBase, sols{1}, 'acetyl-CoA', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'acetate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'CoA', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{2}, 'acetyl-CoA', true, 10^-9);

listMetRxnsWithFluxes(mBase, sols{1}, 'AKG', false, 10^-9);
listMetRxnsWithFluxes(mBase, sols{2}, 'AKG', false, 10^-9);

listMetRxnsWithFluxes(mBase, sols{1}, 'succinate', false, 10^-9);
listMetRxnsWithFluxes(mBase, sols{2}, 'succinate', false, 10^-9);

listMetRxnsWithFluxes(mBase, sols{1}, 'NADH', false, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'NADH', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'ubiquinone', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'ubiquinone', false, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'ATP', false, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'pyruvate', false, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'pyruvate', true, 10^-9);
listMetRxnsWithFluxes(mBase, sols{1}, 'pyruvate', true, 10^-9);

listMetRxnsWithFluxes(mBase, sols{2}, 'succinate', false, 10^-9);


m.ub(strcmp(mBase.rxns, 'MAR09135_REV')) = 1;%L-lactate
%m.ub(strcmp(m.rxns, 'MAR09063_REV')) = 100;%glutamine

%m.ub(strcmp(m.rxns, 'HMR_9133_REV')) = 1000;%pyruvate
%m.ub(strcmp(m.rxns, 'EX_nadh[e]_REV')) = 1000;%NADH









