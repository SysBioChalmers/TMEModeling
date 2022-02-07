%%This file was copied from Gecko and slightly modified.
function OptSigma = modSigmaFitter(model_batch,Ptot,expVal,f)

fprintf('Fitting sigma factor...')

objValues   = [];
errors    = [];
sigParam  = [];
poolIndex = find(strcmpi(model_batch.rxnNames,'prot_pool_exchange'));
objPos    = find(model_batch.c);
%Relax bounds for the objective function
model_batch.lb(objPos) = 0;
model_batch.ub(objPos) = 1000;
lastError = 0;
for i=1:1000
    %Constrains the ecModel with the i-th sigma factor
    sigma = i/1000;
    model_batch.ub(poolIndex) = sigma*Ptot*f;
    solution  = solveLP(model_batch,1);
    if isempty(solution.x)
        solution.x=zeros(length(model_batch.rxns),1);
    end
    objValues = [objValues; solution.x(objPos)];
    error     = ((expVal-solution.x(objPos))/expVal)*100;
    errors    = [errors; error];
    errorText     = num2str(((expVal-solution.x(objPos))/expVal)*100);
    sigParam  = [sigParam; sigma];
    disp([num2str(i) ': ' errorText])
    %stop when we have passed the best value to save computation time
    if (lastError < 0) && (error < 0)
        break;
    end
    lastError = error;
end
fprintf(' Done!\n')
[~, minIndx] = min(errors);
OptSigma     = sigParam(minIndx);
figure
plot(sigParam(1:length(errors)),errors(1:length(errors)),'LineWidth',5)
title('Sigma fitting for growth on glucose minimal media')
xlabel('Average enzyme saturation [-]')
ylabel('Absolute relative error [%]')

end
