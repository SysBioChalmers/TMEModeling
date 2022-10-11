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

rxnsToAdd = struct();
rxnsToAdd.rxns = {'succExp'};
rxnsToAdd.equations = {'succinate[m] => succinate[e]'};
ltModelWithSuccExp = addRxns(ltModel, rxnsToAdd, 3);

tic
succRes = runASimulation(ltModelWithSuccExp, a, bloodData, cell_maintenance);
toc
%save('data/D1_1.mat', 'D1_1')

succRes.resultSolutions{20} %has growth

listMetRxnsWithFluxes(ltModelWithSuccExp, succRes.resultSolutions{20}, 'succinate', false, 10^-9);
%    {'MAR04152'}         0.13731       {'ADP[m] + Pi[m] + succinyl-CoA[m] + 0.0038315 prot_pool[c] => ATP[m] + CoA[m] + succinate[m]'}
%    {'MAR04652'}        0.024165       {'fumarate[m] + ubiquinol[m] + 0.00023282 prot_pool[c] => succinate[m] + ubiquinone[m]'       }
%    {'MAR08743'}       0.0019138       {'FADH2[m] + fumarate[m] + 0.00023282 prot_pool[c] => FAD[m] + succinate[m]'                  }
%    {'succExp' }         0.16339       {'succinate[m] => succinate[e]'   

%We see two things here: 
%1)Reverse complex II is active (small flux)
%2)The TCA cycle up to succinate only can be run, to get more rounds of the TCA cycle per OXPHOS (large flux)

listMetRxnsWithFluxes(ltModelWithSuccExp, succRes.resultSolutions{20}, 'fumarate', false, 10^-9);
%    {'MAR09828'    }      0.00032136       {'fumarate[e] => fumarate[c]'                                           }
%    {'MAR04412_REV'}      0.00011319       {'adenylosuccinate[c] + 0.00044194 prot_pool[c] => AMP[c] + fumarate[c]'}
%    {'MAR04410_REV'}        0.026078       {'malate[m] + 3.6135e-09 prot_pool[c] => fumarate[m] + H2O[m]'          }
%    {'MAR11400_REV'}      0.00032136       {' => fumarate[e]'                                                      }
%so, comes mostly from malate
listMetRxnsWithFluxes(ltModelWithSuccExp, succRes.resultSolutions{20}, 'malate', false, 10^-9);
%    {'MAR04139'    }       0.0001548       {'H+[c] + NADH[c] + OAA[c] + 2.6839e-05 prot_pool[c] => malate[c] + NAD+[c]'}
%    {'MAR04141'    }        0.024336       {'H+[m] + NADH[m] + OAA[m] + 2.6159e-05 prot_pool[c] => malate[m] + NAD+[m]'}
%    {'MAR04408'    }      0.00043184       {'fumarate[c] + H2O[c] + 3.6135e-09 prot_pool[c] => malate[c]'              }
%    {'MAR09030'    }       0.0011641       {'malate[e] => malate[c]'                                                   }
%    {'MAR04852_REV'}       0.0017429       {'AKG[m] + malate[c] + 0.00027666 prot_pool[c] => AKG[c] + malate[m]'       }
%    {'MAR11404_REV'}       0.0011641       {' => malate[e]'                                                            }
%So, mainly from  OAA
listMetRxnsWithFluxes(ltModelWithSuccExp, succRes.resultSolutions{20}, 'OAA', false, 10^-9);
%    {'MAR03827'}       0.0088195       {'AKG[m] + aspartate[m] + 2.9998e-05 prot_pool[c] => glutamate[m] + OAA[m]'                }
%    {'MAR03829'}       0.0001548       {'AKG[c] + aspartate[c] + 3.9896e-05 prot_pool[c] => glutamate[c] + OAA[c]'                }
%    {'MAR04145'}        0.015516       {'citrate[m] + CoA[m] + H+[m] + 8.6015e-05 prot_pool[c] => acetyl-CoA[m] + H2O[m] + OAA[m]'}
listMetRxnsWithFluxes(ltModelWithSuccExp, succRes.resultSolutions{20}, 'citrate', false, 10^-9);
%    {'MAR04458'    }      0.00083364       {'cis-aconitate[m] + H2O[m] + 5.6498e-09 prot_pool[c] => citrate[m]'   }
%    {'MAR05994'    }        0.013773       {'citrate[e] + Na+[e] + 0.00027666 prot_pool[c] => citrate[c] + Na+[c]'}
%    {'MAR04456_REV'}      0.00091328       {'isocitrate[m] + 5.6498e-09 prot_pool[c] => citrate[m]'               }
%    {'MAR09286_REV'}        0.013773       {' => citrate[e]'                                                      }
%    {'MAR02405_REV'}        0.013769       {'AKG[m] + citrate[c] + 0.00027666 prot_pool[c] => AKG[c] + citrate[m]'}


sel = strcmp(ltModelWithSuccExp.rxns,'MAR09415');
sum(sel)
succRes.resultSolutions{20}.x(sel)
%We get succinate export!



