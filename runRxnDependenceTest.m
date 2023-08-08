cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling

%load the "small" model (i.e. just the ltModel, one cell type only)
load('data/ltModel.mat');

cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"

bloodData = prepBloodData();
a = (0.000001:0.000001:0.0001);

res = rxnDependenceTest(ltModel, a(30), a(100), bloodData, cell_maintenance);

%extract the growth rates
hyp = NaN(length(res.hyp),1);
norm = NaN(length(res.norm),1);
for i = 1:length(hyp)
   if ~isempty(res.hyp{i})
       if res.hyp{i}.stat == -1
           hyp(i) = 0;
       else
           hyp(i) = -res.hyp{i}.f;
       end
           
       if res.norm{i}.stat == -1
           norm(i) = 0;
       else
           norm(i) = -res.norm{i}.f;
       end
           
   end
end


sel = ~isnan(hyp);
sum(sel)
ind = find(sel);

hypf = hyp(sel);
normf = norm(sel);

%Now filter all that have a large negative effect for normal cells
normFrac = normf/-res.baseResN.f;
hypFrac = hypf/-res.baseResH.f;

%sel2 = (normFrac > 0.8) & hypFrac < 0.8;
sel2 = normFrac ./ hypFrac > 1.1 & normFrac > 0.5;
sum(sel2)
T = table(ltModel.rxns(ind(sel2)), ...
          normFrac(sel2), ...
          hypFrac(sel2), ...
          constructEquations(ltModel,ltModel.rxns(ind(sel2))))
T.Properties.VariableNames = {'Rxn','NormEff', 'HypoxEff', 'Equation'};
T
writetable(T, 'data/HypoxSensitiveRxns.txt', 'Delimiter', '\t');


%tests
%#1
res.baseResH %-0.0165, fits with the -0.01654112 in Fig. 1, so metabolite setup is ok.
%#2 - Check that all blocked reactions have zero flux, check hyp only, both are done the same way
checked = false(length(res.hyp),1);
checked(cellfun(@isempty, res.hyp)) = true;
for i = 1:length(res.hyp)
    if ~checked(i)
        if res.hyp{i}.stat ~= 1
            checked(i) = true; %infeasible solution, ok
        else
            checked(i) = res.hyp{i}.x(i) == 0;
        end
    end
end
sum(~checked) %0, ok

