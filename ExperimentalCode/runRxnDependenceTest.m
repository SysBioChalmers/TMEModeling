cd C:/Work/MatlabCode/projects/TMEModeling/TMEModeling

%load the "small" model (i.e. just the ltModel, one cell type only)
load('data/ltModel.mat');

cell_maintenance = 1.833; %mmol ATP per gDW and hour, from "A Systematic Evaluation of Methods for Tailoring Genome-Scale Metabolic Models"

bloodData = prepBloodData();
a = (0.000001:0.000001:0.0001);


cd 'ExperimentalCode';
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

sel2 = (normFrac > 0.8) & hypFrac < 0.8;
sum(sel2)
table(ltModel.rxns(ind(sel2)), normFrac(sel2), hypFrac(sel2))

%Let's figure out which combinations we need to test:
sel3 = normFrac < 0.5;
sel3 = hypFrac < 0.5;



