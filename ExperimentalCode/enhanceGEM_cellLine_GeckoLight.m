function enhanceGEM_cellLine_GeckoLight(cellName)
%variant of enhanceGEM_cellLine from the Human-GEM zenodo, adapted to GeckoLight.

%Load original model 
model = load(['../models/' cellName '/' cellName '.mat']);
eval(['model = model.' cellName])
%% model preprocessing 
model_modified = modelModifications(model);
% Save models
save(['../models/' cellName '/model_modified_gl.mat'],'model_modified')
%% Run GeckoLight instead of the GECKO pipeline
%ltModel = CreateECLightModel(model_modified , true, 1);
%run the same way as gecko, i.e., don't set kcat < 1 to 1, and don't fill in missing kcats
ltModel = CreateECLightModel(model_modified , false, 0);

%change the protein pool
% Constrain total protein pool
Ptotal       = 0.593; % Human biomass total protein content [g prot/gDw]
protCoverage = 0.5;
sigma        = 0.5;
constrVal = Ptotal*protCoverage*sigma;

ltModel.ub(length(ltModel.ub)) = constrVal;

% Save output models:
save(['../models/' cellName '/ltModel.mat'],'ltModel')
end
