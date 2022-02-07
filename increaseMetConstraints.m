function exchModel = increaseMetConstraints(model, metaboliteNames, fluxAdditions, checkFeas, includesComp)

if nargin < 5 || isempty(includesComp)
    includesComp = false; %if [s] should be added to the metabolites or not
end

if nargin < 4 || isempty(checkFeas)
    checkFeas = true;
end

[exchRxns, exchIndxs] = getExchangeRxns(model);
exchModel = model;
% exclude protein pool exchange
[~,prot_exch_indx] = ismember({'prot_pool_exchange', 'f_prot_pool_exchange', 'o_prot_pool_exchange'}, model.rxns);
if (length(prot_exch_indx) == 0)
    error('Expected protein pool exchange reaction named "prot_pool_exchange" was not found.');
else
    exchIndxs = setdiff(exchIndxs, prot_exch_indx, 'stable');
    exchRxns  = setdiff(exchRxns, {'prot_pool_exchange', 'f_prot_pool_exchange', 'o_prot_pool_exchange'}, 'stable');
end

% differentiate between uptake and production reactions
uptkIndxs = exchIndxs(contains(exchRxns, '_REV'));

% open uptake of media components one by one
unusedMets = [];
for i = 1:length(metaboliteNames)

    % get metabolite index
    if (includesComp)
        metIndx = getIndexes(model, metaboliteNames{i}, 'metcomps');
    else
        metIndx = getIndexes(model, strcat(metaboliteNames{i},'[s]'), 'metcomps');
    end

    % get rxns for metabolite
    metRxns = find(model.S(metIndx,:));

    % get the uptake reaction for the metabolite
    metUptakeRxn = intersect(metRxns,uptkIndxs);
    if isempty(metUptakeRxn)
        unusedMets = [unusedMets; metaboliteNames(i)];
    else
        exchModel.ub(metUptakeRxn) = exchModel.ub(metUptakeRxn) + fluxAdditions(i); %important to use exchModel on both sides here (and not model), since metabolites can appear more than once in the list
    end
end


% report unused metabolites
if ~isempty(unusedMets)
    fprintf('WARNING: The following metabolites are either not in the model or do not have exchange reactions:\n');
    fprintf('\t%s\n',unusedMets{:});
end

%Check if model is feasible
if ( checkFeas )
    sol = solveLP(exchModel);
    if ~isempty(sol.x)
        disp(['Constrained ec model "' exchModel.id '" is feasible'])
    else
        disp(['*** Constrained ec model "' exchModel.id '" is INFEASIBLE ***'])
        exchModel = [];
    end
end
