function outModel = fixFluxes(inModel, solution, rxns, setUpperBound, margin)
% 
% Sets the upper or lower bound of reactions to the value given by the
% solution (with some minor margin). Typically, you should first solve with
% these reactions in the objective function, then use this function to fix
% them, and then optimize on something else, all to avoid randomness in the final
% solution.
% 
%
% Input:
%
%
% Output:
%
%

if nargin < 5
    margin = 10^3;
end

outModel = inModel;


for i = 1:length(rxns)
   marginVal = solution.x(strcmp(outModel.rxns,rxns{i})) / margin;
   if marginVal < 10^-10
       marginVal = 10^-10;
   end
   if (setUpperBound)
      outModel.ub(strcmp(outModel.rxns,rxns{i})) = solution.x(strcmp(outModel.rxns,rxns{i})) + marginVal;
   else
      outModel.lb(strcmp(outModel.rxns,rxns{i})) = solution.x(strcmp(outModel.rxns,rxns{i})) - marginVal;
   end
end