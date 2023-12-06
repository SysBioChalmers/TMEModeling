%This is a modified version of predict_cellLines adapted for GeckoLight
cd C:/Work/MatlabCode/projects/HumanGemZenodo/Human1_Publication_Data_Scripts/ec_GEMs/ComplementaryScripts

cellNames   = [{'HS_578T'} {'RPMI_8226'} {'HT29'} {'MALME_3M'} {'SR'}...
               {'UO_31'} {'MDMAMB_231'} {'HOP62'} {'NCI_H226'} {'HOP92'} {'O_786'}];
constraints = [{'media'} {'glucose'} {'L-lactate'} {'threonine'}];
         
pred_Grates    = [];
pred_EC_Grates = [];
meanError      = [];
mean_EC_error  = [];
for constLevel = [0 1 2 3]
    %Run growth simulations for ecModels
    [meanErr,pred,meas_GRates] = ExchFluxesComparison_NCI60_GeckoLight(constLevel,true,false);
    pred_EC_Grates   = [pred_EC_Grates, pred];
    mean_EC_error    = [mean_EC_error; meanErr]; 
    %Run growth simulations for standard models
    [meanErr,pred,~] = ExchFluxesComparison_NCI60_GeckoLight(constLevel,false,false);
    pred_Grates      = [pred_Grates, pred];
    meanError        = [meanError; meanErr]; 
end
experimental = meas_GRates';
cellNames    = cellNames';
%Write results ecModels
results    = table(cellNames,experimental);
errorTable = table(cellNames);
for i=1:length(constraints)
    results    = [results table(pred_EC_Grates(:,i),'VariableNames',constraints(i))];
    errorVec   = abs(experimental-pred_EC_Grates(:,i))./experimental;
    errorTable = [errorTable table(errorVec,'VariableNames',constraints(i))];
end
fileName = '../Results/11_cellLines_NCI60/ecModels_gRates.txt';
writetable(results,fileName,'Delimiter','\t','QuoteStrings',false);
fileName = '../Results/11_cellLines_NCI60/ecModels_error_gRates.txt';
writetable(errorTable,fileName,'Delimiter','\t','QuoteStrings',false);
%Write results Models
results  = table(cellNames,experimental);
errorTable = table(cellNames);
for i=1:length(constraints)
    results    = [results table(pred_Grates(:,i),'VariableNames',constraints(i))];
    errorVec   = abs(experimental-pred_Grates(:,i))./experimental;
    errorTable = [errorTable table(errorVec,'VariableNames',constraints(i))];
end
fileName = '../Results/11_cellLines_NCI60/models_gRates.txt';
writetable(results,fileName,'Delimiter','\t','QuoteStrings',false);
fileName = '../Results/11_cellLines_NCI60/models_error_gRates.txt';
writetable(errorTable,fileName,'Delimiter','\t','QuoteStrings',false);
