%This code was used to estimate a reasonable protein usage constraint value in GECKO Light
%by fitting the value to NCI 60 growth data.
cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling
load('data/ltModel.mat')

%read experimental data
expData = readtable('data/exchangeFluxes_NCI60_Conv.txt');

ptotAll = table2array(expData(strcmp(expData.metabolite, 'Ptot'),3:end));
growthAll = table2array(expData(strcmp(expData.metabolite, 'biomass'),3:end));
nn = length(ptotAll);
sigmas = NaN(nn,1);

metInd = 1:length(expData.RxnID) - 2;

tableVals = table2array(expData(:,3:end));


for i = 1:length(ptotAll)
    disp(['Running cell line: ' num2str(i)])
    %model = prepMinModel;
    model = ltModel;
    model = setGrowthMedium(model, true, 'Hams');
    model.lb(model.lb == -1000) = -Inf;
    model.ub(model.ub == 1000) = Inf;

    constr = tableVals(metInd,i);

    %It seems that some essential amino acids are sometimes limiting. This can be
    %due to for example varying amino acid composition in the biomass. Relax that
    %constraint.
    %Set values lower than 0.03 to 0.03 for essential amino acids
    rge = [9 10 12 14 16 21 23 25];%histidine not included in table
    %expData(rge,2)
    rgeChange = rge(constr(rge) > -0.03);
    constr(rgeChange) = -0.03;

    sel = constr < 0;
    constr(sel) = -constr(sel);
    rxnNames = expData.RxnID(metInd);
    rxnNames(sel) = strcat(rxnNames(sel), '_REV');
    
    rxnsToZero = expData.RxnID(metInd);
    rxnsToZero(~sel) = strcat(rxnsToZero(~sel), '_REV');
    
    %let's just loop to be safe
    for j=1:length(rxnNames)
        %Relax the constraints to 1.5 of the measured values. This is to
        %account for differences in kcats, biomass reaction etc between
        %the real cells and the model.
        %The constraints are not that important, the cells are taking up as
        %much as they can - the most important thing is to not constrain it
        %too hard
        %or, just don't...
        model.ub(strcmp(model.rxns, rxnNames(j))) = constr(j);%*1.5; 
    end
    for j=1:length(rxnsToZero)
        %we have problems with folate if we constrain the reverse reactions, so
        %skip that
        if (~strcmp(rxnsToZero(j), 'MAR09146_REV'))
            model.ub(strcmp(model.rxns, rxnsToZero(j))) = 0;
        end
    end
    lim = 0.06;
    model.ub(strcmp(model.rxns, 'MAR09423_REV')) = lim; %thymidine
    
    model.c = double(strcmp(model.rxns,'MAR13082'));%set objective function to biomass_human
    %add NGAM
    model.lb(strcmp(model.rxns,'MAR03964')) = 1.833; %non-growth associated cell maintenance
    
    sigmas(i) = modSigmaFitter(model,ptotAll(i),growthAll(i),0.5); %we just use f = 0.5
    %solution  = solveLP(model,1);

end


proBounds = NaN(length(sigmas),1);

for i=1:length(sigmas)
    proBounds(i) = sigmas(i)*ptotAll(i)*0.5;
end
proBounds
figure
scatter(growthAll, proBounds)

mean(proBounds)%0.0224 - this is the value used in GECKO Light
sprintf('%d',mean(proBounds))%2.238315e-02
sprintf('%d',max(proBounds))%2.798355e-02

save('data/FitProteinConstraintResults.mat', 'sigmas', 'ptotAll', 'proBounds', 'growthAll')
%load('data/FitProteinConstraintResults.mat')


