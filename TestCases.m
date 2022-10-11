cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling %the folder with the functions



%% TC002 - buildFullModel2 + addCollaborationCost + blockCollaboration
load('data/ltModel.mat');



%Generate a combination of runs with fracECM = 0.01, 0.25, 0.5 - fracC = 0.6, fracF = 0.2
%Then vary frac fibroblast = 0.05, 0.2, 0.4 - fracC is adapted, fracECM = 0.25
%Note that the middle in both are the same, so we get 5 runs in total
%Run each with both type of runs
mx = addCollaborationCost(buildFullModel2(ltModel, 0.6, 0.2, 0.01, 0.2));

%check the total_biomass function
constructEquations(mx, 'total_biomass') %{'biomass[c] + 0.016835 ECM_biomass[f_e] => '}
%the stochiomentry looks reasonable
constructEquations(mx, 'ECM_biomass') %{'0.063828 heparan sulfate proteoglycan[f_e] + 8.4379 ECM_protein_pool_biomass[f_e] => ECM_biomass[f_e]'}
%looks reasonable - no point in checking stochiometry, that is calculated and just entered
constructEquations(mx, 'ecm_protein_pool_biomass') %{'0.271 glycyl-tRNA(gly)[f_c] + 0.0948 L-alanyl-tRNA(ala)[f_c] + 0.0499 L-arginyl-tRNA(arg)[f_c] + 0.0228 L-asparaginyl-tRNA(asn)[f_c] + 0.0405 L-aspartyl-tRNA(asp)[f_c] + 0.0104 L-cysteinyl-tRNA(cys)[f_c] + 0.0304 L-glutaminyl-tRNA(gln)[f_c] + 0.0503 L-glutamyl-tRNA(glu)[f_c] + 0.0078 L-histidyl-tRNA(his)[f_c] + 0.0187 L-isoleucyl-tRNA(ile)[f_c] + 0.0367 L-leucyl-tRNA(leu)[f_c] + 0.0382 L-lysyl-tRNA(lys)[f_c] + 0.0084 L-methionyl-tRNA(met)[f_c] + 0.0177 L-phenylalanyl-tRNA(phe)[f_c] + 0.1832 L-prolyl-tRNA(pro)[f_c] + 0.04 L-seryl-tRNA(ser)[f_c] + 0.0307 L-threonyl-tRNA(thr)[f_c] + 0.004 L-tryptophanyl-tRNA(trp)[f_c] + 0.0098 L-tyrosyl-tRNA(tyr)[f_c] + 0.0348 L-valyl-tRNA(val)[f_c] => 0.0948 tRNA(ala)[f_c] + 0.0499 tRNA(arg)[f_c] + 0.0228 tRNA(asn)[f_c] + 0.0405 tRNA(asp)[f_c] + 0.0104 tRNA(cys)[f_c] + 0.0304 tRNA(gln)[f_c] + 0.0503 tRNA(glu)[f_c] + 0.271 tRNA(gly)[f_c] + 0.0078 tRNA(his)[f_c] + 0.0187 tRNA(ile)[f_c] + 0.0367 tRNA(leu)[f_c] + 0.0382 tRNA(lys)[f_c] + 0.0084 tRNA(met)[f_c] + 0.0177 tRNA(phe)[f_c] + 0.1832 tRNA(pro)[f_c] + 0.04 tRNA(ser)[f_c] + 0.0307 tRNA(thr)[f_c] + 0.004 tRNA(trp)[f_c] + 0.0098 tRNA(tyr)[f_c] + 0.0348 tRNA(val)[f_c] + ECM_protein_pool_biomass[f_c]'}
%looks reasonable - no point in checking stochiometry, that is calculated and just entered

%check that communication from other cells is blocked - use glucose as example
scomp = find(strcmp(mx.comps,'e'));
SSub = mx.S(strcmp(mx.metNames,'glucose') & (mx.metComps == 1),:);
size(SSub) %looks good, just one row
%now list the reactions that are non-zero in SSub
rxnsOfInterestSel = (SSub ~= 0).';
constructEquations(mx, mx.rxns(rxnsOfInterestSel) )
%removed irrelevant lines
%    {'glucose[c] => glucose[e]'                                                                                           }
%    {'glucose[e] => '                                                                                                     }
%    {'glucose[e] => glucose[c]'                                                                                           }
%    {' => glucose[e]'                                                                                                     }
%    {'1e-06 prot_pool[c] + glucose[f_e] => glucose[e]'                                                                    } - transport from fibroblasts, ok
%    {'glucose[e] => glucose[f_e]'                                                                                         }
%    {'glucose[e] => glucose[o_e]' 

%Check: there is communication to f_s and o_s, and back from f_s, but at a cost (for addCollaborationCost), but not back from o_s
%so, all is as it should be.

%now check blockCollaboration - in the example below, '1e-06 prot_pool[c] + glucose[f_e] => glucose[e]' should change to '1e-06 prot_pool[c] + glucose[f_e] => '
mxb = blockCollaboration(mx, 'L-lactate');
constructEquations(mxb, mxb.rxns(rxnsOfInterestSel))
%reaction changes to this, ok
%{'1e-06 prot_pool[c] + glucose[f_e] => '                                                                             }




%% TC003 - 
load('data/ltModel.mat');
modelTmp = constrainMedium(ltModel, {'glucose[e]'}, [7 13], false, true);

indices = find(modelTmp.ub ~= 0)
[exchRxns, exchIndxs] = getExchangeRxns(modelTmp,'in');

exchIndxsFilt = exchIndxs(modelTmp.ub(exchIndxs) ~= 0);
modelTmp.ub(exchIndxsFilt) %glucose changed to 13, ok
constructEquations(modelTmp, modelTmp.rxns(exchIndxsFilt)) %in total 45, all fully open except prot_pool
%so, setGrowthMedium has implicitly been tested here as well, since it is run on ltModel

%now check that we don't leave open metabolites that shouldn't be open
bloodData = prepBloodData();
ux = zeros(length(bloodData.totDxC),2);
ux(:,2) = bloodData.totDxC*0.0001;
%insert constraints from blood conc + Ham's media
modelTmp2 = constrainMedium(ltModel, bloodData.totMets, ux, false, true);
indices = find(modelTmp2.ub ~= 0);
[exchRxns, exchIndxs] = getExchangeRxns(modelTmp2,'in');
exchIndxsFilt = exchIndxs(modelTmp2.ub(exchIndxs) ~= 0);

constructEquations(modelTmp2, modelTmp2.rxns(exchIndxsFilt)) %in total 91
%check which ones are fully open
exchIndxsFilt2 = exchIndxs(modelTmp2.ub(exchIndxs) >= 1000);
constructEquations(modelTmp2, modelTmp2.rxns(exchIndxsFilt2)) 
%    {' => H2O[s]'    }
%    {' => Pi[s]'     }
%    {' => sulfate[s]'}
%    {' => Fe2+[s]'   }
% OK


%TC004 - HaveFluxOptimized
%Run both HaveFlux and HaveFluxOptimized and see that they produce identical results.
%Do so for a smaller model, just to keep it smaller
load('data/ltModel.mat');

bloodData = prepBloodData();
ux = zeros(length(bloodData.totDxC),2);
ux(:,2) = bloodData.totDxC*0.0001;
%insert constraints from blood conc + Ham's media
modelTmp = constrainMedium(ltModel, bloodData.totMets, ux, false, true);
x = haveFluxOptimized(modelTmp);
y = haveFlux(modelTmp);
all(x == y)%false - something is not right perhaps?
find(x ~= y)
x(x ~= y)
y(x ~= y)
%There are a few reactions that differ - investigate
rxnsThatDiffer = modelTmp.rxns(x ~= y);
%debugging haveFlux
haveFlux(modelTmp, 10^-14, rxnsThatDiffer)
haveFluxOptimized(modelTmp, 10^-14, rxnsThatDiffer)


%TC005 - FixFluxes
load('data/ltModel.mat');
bloodData = prepBloodData();
ux = zeros(length(bloodData.totDxC),2);
ux(:,2) = bloodData.totDxC*0.0001;
%insert constraints from blood conc + Ham's media
modelTmp2 = constrainMedium(ltModel, bloodData.totMets, ux, false, true);
modelTmp2.c = double(strcmp(modelTmp2.rxns, 'MAR13082'));
res = solveLP(modelTmp2,1);
modelTmp3 = fixFluxes(modelTmp2, res, {'MAR13082'}, false);
(-res.f) - modelTmp3.lb(strcmp(modelTmp3.rxns,'MAR13082'))%8.8359e-05
(-res.f)%0.0883
%So, the margin is much smaller, ok.

%TC006 - runMacrophageSimulation + increaseMetConstraints
%First, remove comments from the test code in the runMacrophageSimulation function
%Then set a breakpoint there, then run
load('data/ltModelFull.mat');
macrData = readtable('data/MacrophageInput.txt');
bloodData = prepBloodData();
tmp = runMacrophageSimulation(ltModelFull, 0.0002, bloodData, macrData, 0.1); %important to run on the full model - there are new exchange reactions used here that were not used when minimizing the model
%check that the output from expRes - tmp is really small, should only be a small roundoff error

%TC007 - listMetRxnsWithFluxes
load('data/ltModel.mat');
bloodData = prepBloodData();
ux = zeros(length(bloodData.totDxC),2);
ux(:,2) = bloodData.totDxC*0.0001;
%insert constraints from blood conc + Ham's media
modelTmp2 = constrainMedium(ltModel, bloodData.totMets, ux, false, true);
modelTmp2.c = double(strcmp(modelTmp2.rxns, 'MAR13082'));
res = solveLP(modelTmp2,1);
listMetRxnsWithFluxes(modelTmp2, res, 'glucose', true, 10^-9);
%    {'MAR07747'    }        3.5643         {'ADP[c] + glucose[c] + 1.2711e-05 prot_pool[c] => AMP[c] + glucose-6-phosphate[c] + H+[c]'}
%    {'MAR05029_REV'}        3.5643         {'glucose[e] => glucose[c]'                                                                }
%as expected, glucose is taken up and enters glycolysis, all reactions have glucose on the left side
listMetRxnsWithFluxes(modelTmp2, res, 'glucose', false, 10^-9);
%    {'MAR05029_REV'}        3.5643         {'glucose[e] => glucose[c]'}
%    {'MAR09034_REV'}        3.5643         {' => glucose[e]'          }
% looks reasonable, here we get the glucose reactions with glucose on the right side, different compartments covered
listMetRxnsWithFluxes(modelTmp2, res, 'glucose', false, 4); %here they are filtered, ok
listMetRxnsWithFluxes(modelTmp2, res, 'glucose', false, 2); %here they are not filtered, ok
%So ok, the function seems to work

%also test with compartment
listMetRxnsWithFluxes(modelTmp2, res, 'L-lactate[e]', true, 10^-9);
%{'MAR09135'}        7.0234         {'L-lactate[e] => '} 
%ok
