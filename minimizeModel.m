function outModel = minimizeModel(ecModelCorr, bloodData)
% 
% Removes unnecessary reactions from the model, such as dead end rxns etc. Takes some time to run (> 1 hour).
% 

outModel = ecModelCorr;

outModel = setGrowthMedium(outModel, true, 'Hams');
%need to do this afterwards
outModel.lb(outModel.lb == -1000) = -Inf;
outModel.ub(outModel.ub == 1000) = Inf;

ux = zeros(length(bloodData.totDxC),2);
a = 0.0002;
ux(:,2) = bloodData.totDxC*a;
%insert constraints from blood conc + Ham's media
outModel = constrainMedium(outModel, bloodData.totMets, ux, false, true);


hfRes = haveFlux(outModel,10^-14);
save('data/hfRes.mat', 'hfRes'); %just to be safe, takes ~24 hours to run

%before: S: [15416×46309 double]
outModel = removeReactions(outModel,~hfRes, true, true, false); %don't remove compartments, could perhaps lead to problems later.
%after: S: [13500×39240 double]
