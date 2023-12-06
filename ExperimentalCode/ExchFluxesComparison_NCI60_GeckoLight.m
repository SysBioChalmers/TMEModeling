function [meanGRerror,pred_GRates,meas_GRates] = ExchFluxesComparison_NCI60_GeckoLight(constLevel,ecFlag,fixedBounds)
%Variant of ExchFluxesComparison_NCI60 adapted for GeckoLight
clc
close all
current      = pwd;
%Assign specific color code for the different cell-lines
[cellLines, colorS] = assignNamesAndColors;
%Open exchange fluxes data file
fID     = fopen('../DataFiles/exchangeFluxes_NCI60.txt');
expData = textscan(fID,'%s %s %f %f %f %f %f %f %f %f %f %f %f',...
                   'Delimiter','\t','HeaderLines',1);
%%Save experimental data 
expIDs = expData{1};expMets = expData{2};expData = expData(3:end);
%Initialize variables
legendStr       = {};
errors_ExFlux   = zeros(length(cellLines),3);
errors_GRates   = zeros(length(cellLines),1);
pred_GRates     = zeros(length(cellLines),1);
errors_ExFluxLt   = zeros(length(cellLines),3);
errors_GRatesLt   = zeros(length(cellLines),1);
pred_GRatesLt     = zeros(length(cellLines),1);
%Create output files name strings
if ecFlag
    fileName1 = ['ecModels_const_' num2str(constLevel) '_exchangeFluxesComp.txt'];
    fileName2 = ['ecModels_const_' num2str(constLevel) '_errorMetrics.txt'];
    fileName3 = ['ecModels_const_' num2str(constLevel) '_error_GRate.txt'];
else
    fileName1 = ['models_const_' num2str(constLevel) '_exchangeFluxesComp.txt'];
    fileName2 = ['models_const_' num2str(constLevel) '_errorMetrics.txt'];
    fileName3 = ['models_const_' num2str(constLevel) '_error_GRate.txt'];
end

%For each cell line
for i=1:length(cellLines)
    %Measured fluxes for the i-th cell-line
    measuredFluxes = expData{i};
    cellLineStr    = strrep(cellLines{i},'_','-');
    mkdir (['../models/' cellLines{i} '/Results'])
    disp(cellLines{i})
    %Load models
    load(['../models/' cellLines{i} '/ecModel_batch.mat'])
    load(['../models/' cellLines{i} '/ltModel.mat'])
    load(['../models/' cellLines{i} '/' cellLines{i} '.mat'])
    eval(['model =' cellLines{i} ';'])
    model.b = zeros(length(model.mets),1);
    %Allow any level of average enzyme saturaion (0 - 100%)
    P_index = strcmpi(ecModel_batch.rxnNames,'prot_pool_exchange');
    ecModel_batch.ub(P_index) = 0.593;
    %Save models
    %save (['../models/' cellLines{i} '/ecModel_batch.mat'], 'ecModel_batch')
    %save (['../models/' cellLines{i} '/model.mat'], 'model')
    %Run FBA 
    ecMsol = solveLP(ecModel_batch);
    ltMsol = solveLP(ltModel);
    GEMsol = solveLP(model);
    %If both models are feasible
    if ~isempty(ecMsol.x) & ~isempty(ltMsol.x) & ~isempty(GEMsol.x)
        cd Simulation
        %Set media composition
        ecModel_batch  = setHamsMedium_GeckoLight(ecModel_batch,true);
        ltModel  = setHamsMedium_GeckoLight_gl(ltModel,true);
        model          = setHamsMedium_GeckoLight(model,false);
        %Set fixed constraints (nutrients uptakes) 
        ecModel_batch = setDataConstraints(ecModel_batch,expData{i},expIDs,true,constLevel,fixedBounds,false);
        ltModel = setDataConstraints(ltModel,expData{i},expIDs,true,constLevel,fixedBounds,true);
        model         = setDataConstraints(model,expData{i},expIDs,false,constLevel,fixedBounds,false);
        %Get exchange fluxes and metabolites IDs in the original model
        %model = rmfield(model,'unconstrained');
        [EX_IDs,exc_Indexes] = getExchangeRxns(model);
        exchangeMets = {};
        for j=1:length(expIDs)-1
            metJ         = strcmpi(model.rxns,expIDs(j));
            exchangeMets = [exchangeMets;...
                model.metNames(find(model.S(:,metJ),1))];
        end
        %Get parsimonious solutions for both models
        ecMsol = minProtSimulation_GeckoLight(ecModel_batch);
        ltMsol = minProtSimulation_GeckoLight(ltModel);
        GEMsol = solveLP(model,1);
        %If both models are feasible then save results
        if ~isempty(ecMsol.x) & ~isempty(ltMsol.x) & ~isempty(GEMsol.x)
            %Extract values for exchange fluxes from GEMsol
            exc_Model = GEMsol.x(exc_Indexes);
            cd ..
            exc_ecModel =[];
            ubEC = [];
            lbEC = [];
            %Map each exchange rxn to ecModel_batch
            for j=1:length(exc_Indexes)
                indexes = rxnMapping(EX_IDs{j},ecModel_batch,true);
                if length(indexes)>1
                    exc_ecModel = [exc_ecModel; ecMsol.x(indexes(1))-ecMsol.x(indexes(2))];
                    ubEC = [ubEC; ecModel_batch.ub(indexes(1))+ecModel_batch.ub(indexes(2))];
                    lbEC = [lbEC; ecModel_batch.lb(indexes(1))+ecModel_batch.lb(indexes(2))];
                else
                    exc_ecModel = [exc_ecModel; ecMsol.x(indexes)];
                    ubEC = [ubEC; ecModel_batch.ub(indexes(1))];
                    lbEC = [lbEC; ecModel_batch.lb(indexes(1))];
                end
            end
            exc_ltModel =[];
            ubLt = [];
            lbLt = [];
            %Map each exchange rxn to ltModel
            for j=1:length(exc_Indexes)
                indexes = rxnMapping(EX_IDs{j},ltModel,true);
                if length(indexes)>1
                    exc_ltModel = [exc_ltModel; ltMsol.x(indexes(1))-ltMsol.x(indexes(2))];
                    ubLt = [ubLt; ltModel.ub(indexes(1))+ltModel.ub(indexes(2))];
                    lbLt = [lbLt; ltModel.lb(indexes(1))+ltModel.lb(indexes(2))];
                else
                    ubLt = [ubLt; ecModel_batch.ub(indexes(1))];
                    lbLt = [lbLt; ecModel_batch.lb(indexes(1))];
                    exc_ltModel = [exc_ltModel; ltMsol.x(indexes)];
                end
            end
            
            %Keep just the exchange fluxes that are also part of the
            %dataset and compare both GEM and ecGEM predictions
            [~,toKeep,order]  = intersect(EX_IDs,expIDs);
            if length(toKeep) == length(expIDs(1:end-1))
                if ecFlag
                    predictions = exc_ecModel(toKeep);
                    predictionsLt = exc_ltModel(toKeep);
                    fileName1 = ['ecModels_const_' num2str(constLevel) '_exchangeFluxesComp.txt'];
                    fileName2 = ['ecModels_const_' num2str(constLevel) '_errorMetrics.txt'];
                    fileName3 = ['ecModels_const_' num2str(constLevel) '_error_GRate.txt'];
                    fileName1Lt = ['ecModels_const_' num2str(constLevel) '_exchangeFluxesComp_Lt.txt'];
                    fileName2Lt = ['ecModels_const_' num2str(constLevel) '_errorMetrics_Lt.txt'];
                    fileName3Lt = ['ecModels_const_' num2str(constLevel) '_error_GRate_Lt.txt'];
                else
                    predictions = exc_Model(toKeep);
                    predictionsLt = exc_ltModel(toKeep);
                    fileName1 = ['models_const_' num2str(constLevel) '_exchangeFluxesComp.txt'];
                    fileName2 = ['models_const_' num2str(constLevel) '_errorMetrics.txt'];
                    fileName3 = ['models_const_' num2str(constLevel) '_error_GRate.txt'];
                end
                experimental  = measuredFluxes(order);
                exchangeMets  = exchangeMets(order);
                exchangeIDs   = expIDs(order);
                if i==1
                    predictedFluxes = table(exchangeIDs,exchangeMets);
                end
                varStr = [{['exp_' cellLines{i}]} cellLines(i), {['lt_' cellLines{i}]}];
                tempT  = table(experimental,predictions,predictionsLt,'VariableNames',varStr);
                predictedFluxes = [predictedFluxes tempT];
                %Get biomass exchange index
                bioIndx = find(strcmpi(exchangeMets,'biomass'));
                [Xvalues,Yvalues,direction,errors_ExFlux(i,:)] = getPlotValues(experimental,predictions);
                %Calculate mean absolute error for Growth rate predictions
                errors_GRates(i) = computeErrorMetric(experimental(bioIndx),predictions(bioIndx),'MRE');
                pred_GRates(i)   = predictions(bioIndx);
                meas_GRates(i)   = experimental(bioIndx);
                disp(['Relative error for GRate prediction: ' num2str(errors_GRates(i))])
                plotDataPoints(Xvalues,Yvalues,colorS(i,:),direction)
                legendStr = [legendStr; [cellLineStr '/ RSME = ' num2str(errors_ExFlux(i,3))]];
                hold on
                
                %Calculate mean absolute error for Growth rate predictions - Lt
                [Xvalues,Yvalues,direction,errors_ExFluxLt(i,:)] = getPlotValues(experimental,predictionsLt);
                errors_GRatesLt(i) = computeErrorMetric(experimental(bioIndx),predictionsLt(bioIndx),'MRE');
                pred_GRatesLt(i)   = predictions(bioIndx);
            else
                disp('Inconsistent mapping')
            end
        else
            disp(['Not feasible 2: ' cellLines{i}])
        end
    else
        disp(['Not feasible 1: ' cellLines{i}])
    end
end

%Plot exchange fluxes prediction comparison
x1 = linspace(-9,1,100);
plot(x1,x1)
legend(legendStr)
hold on
%Write output files
mkdir ../Results/11_cellLines_NCI60
cd ../Results/11_cellLines_NCI60

%%write file with different error metrics for all cell lines
%variables = {'cell_Line' 'pearson' 'MAE' 'RMSE'};
%T         = table(cellLines',errors_ExFlux(:,1),errors_ExFlux(:,2),errors_ExFlux(:,3),'VariableNames',variables);
%writetable(T,fileName2,'QuoteStrings',false,'Delimiter','\t')

%%write file with GRate predictions MAE
%variables = {'cell_Line' 'MAE'};
%T         = table(cellLines',errors_GRates,'VariableNames',variables);
%writetable(T,fileName3,'QuoteStrings',false,'Delimiter','\t')

%%write file with compared exchange fluxes for all cell-lines
writetable(predictedFluxes,fileName1,'QuoteStrings',false,'Delimiter','\t')
meanGRerror = mean(errors_GRates);
% minGRerror  = min(errors_GRates);
% maxGRerror  = max(errors_GRates);
cd (current)
end
%--------------------------------------------------------------------------
function [cellNames, colors] = assignNamesAndColors
cellNames = [{'HS_578T'} {'RPMI_8226'} {'HT29'} {'MALME_3M'} {'SR'}...
           {'UO_31'} {'MDMAMB_231'} {'HOP62'} {'NCI_H226'} {'HOP92'}...
           {'O_786'}];

colors    = [0    0    128             %dark blue
             0    0    255             %blue
             0    170  255             %light blue
             0    128  0               %forest green
             128  255  0               %lime green
             255  105  180             %pink
             255  200  050             %yellow
             191  0    255             %purple
             255  128  0               %orange
             0    0    0               %black
             255  0    0]./255;        %red
end     
%--------------------------------------------------------------------------
function [Xvalues,Yvalues,direction,errors] = getPlotValues(predictions,experimental)
Xvalues      = (log10(abs(experimental)+1E-9));
Yvalues      = (log10(abs(predictions)+1E-9));
direction    = sign(experimental);
MRE          = computeErrorMetric(experimental,predictions,'MRE');
RMSE         = computeErrorMetric(experimental,predictions,'RMSE');
pearson      = computeErrorMetric(experimental,predictions,'pearson');
errors       = [pearson,MRE,RMSE];
end
%--------------------------------------------------------------------------
function plotDataPoints(Xvalues,Yvalues,color,direction)
sizes   = [];
colorS  = zeros(length(Xvalues),3);
for i=1:length(Xvalues)
    if direction(i)<=0
        sizes = [sizes; 75];
    else
        sizes = [sizes; 150];
    end
    colorS(i,:)  = color;
end
scatter(Xvalues,Yvalues,sizes,colorS,'d','filled')
end
%--------------------------------------------------------------------------
function avg_Eps = computeErrorMetric(data,predictions,method)
data(data==0) = 1E-9;
switch method
    case 'pearson'
        avg_Eps  = corrcoef(predictions,data);
        avg_Eps  = avg_Eps(1,2);
    case 'MRE'
        relErr   = (data-predictions)./data;
        relErr   = (abs(relErr));
        avg_Eps  = mean(relErr);
    case 'RMSE'
        total = 0;
        T     = length(data);
        for i=1:T
            total = total + (data(i) - predictions(i))^2;
        end
        avg_Eps = sqrt(total/(T));
end
end
%--------------------------------------------------------------------------
function  solution = minProtSimulation_GeckoLight(model)
sol      = solveLP(model,1);
objIndex = find(model.c);
objValue = sol.x(objIndex);
%Fix objective value
model.ub(objIndex) = 1.001*objValue;
model.lb(objIndex) = 0.999*objValue;
%Find protein pool exchange pseudoreaction and set it as an objetive to
%minimize
model.c           = zeros(length(model.c),1);
protIndx          = find(strcmpi(model.rxns,'prot_pool_exchange'));
model.c(protIndx) = -1;
solution          = solveLP(model,1);
end
%--------------------------------------------------------------------------
function model = setDataConstraints(model,fluxes,expIDs,ecFlag,const_level,fixed,light)
%Biomass UB
GrowthRate = fluxes(end-1);
%value   = GrowthRate;
GrowthIndx = find(strcmpi(model.rxns,expIDs(end-1)));
%model   = setBounds(model,GrowthIndx,value,ecFlag,true);
if ecFlag
    Prot_biomass = fluxes(end); %Total protein content in biomass [g prot/g DW]
    Prot_pool    = fluxes(end); %Amount of enzymes available for biochemical reactions
    model        = rescaleBiomassProtein(model,Prot_pool,ecFlag,light);
end
%Glucose exchange bound
if const_level >0
    %value   = GrowthRate;
    %model   = setBounds(model,GrowthIndx,value,ecFlag,true);
    index   = strcmpi(expIDs,'HMR_9034');
    rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
    value   = fluxes(index);
    model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
    %L-lactate exchange bound
    if const_level>1
        index   = strcmpi(expIDs,'HMR_9135');
        rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
        value   = fluxes(index);
        model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
        %Threonine exchange bound
        if const_level>2
            index   = strcmpi(expIDs,'HMR_9044');
            rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
            value   = fluxes(index);
            model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
        end
    end
end
end
%--------------------------------------------------------------------------
function model = setBounds(model,index,value,ecFlag,fixed)
if fixed
    lowValue = 0.9;
else
    lowValue = 0;
end
direction = sign(value);
if direction <0
    if ecFlag
        %Block the opposite direction (production)
        model.ub(index) = 0;
        %Find uptake reaction
        rxn   = [model.rxns{index} '_REV'];
        index = find(strcmpi(model.rxns,rxn));
        model.ub(index) = abs(value);
        model.lb(index) = lowValue*abs(value);
    else
        model.lb(index) = value;
        model.ub(index) = lowValue*value;
    end
else
    model.ub(index) = value;
    model.lb(index) = lowValue*value;
end
end
%--------------------------------------------------------------------------
function model = rescaleBiomassProtein(model,Ptot,constrainPool,light)
if nargin <3 
    gIndex    = find(strcmpi(model.rxns,'biomass_human'));
    protIndex = find(strcmpi(model.metNames,'protein_pool_biomass'));
    coeff     = model.S(protIndex,gIndex); %Extract coeff corresponding to 0.593 g Prot / g Biomass
    newCoeff  = coeff*Ptot/0.593;
    %Rescale stoichiometric coefficient of proteinPool in biomass equation
    model.S(protIndex,gIndex) = newCoeff;
else
    if constrainPool
        protPool           = find(strcmpi(model.rxns,'prot_pool_exchange'));
        if light
            model.ub(protPool) = Ptot*0.5;
        else
            model.ub(protPool) = Ptot;
        end
    end
end
end
