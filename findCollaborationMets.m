%m2, a, bloodData, cell_maintenance
%inModel = m2;
%cellMaintenance = cell_maintenance;
function res = findCollaborationMets(inModel, a, bloodData, cellMaintenance)
nPoints = length(a);
res.succeeded = false(nPoints,1);
res.a = a;
res.collaborationMets = false(length(inModel.metNames),nPoints);
res.rxns = inModel.rxns;
res.metNames = inModel.metNames;

for i = 1:nPoints
    disp(['A index: ', num2str(i)])
    ux = zeros(length(bloodData.totDxC),2);
    ux(:,2) = bloodData.totDxC*a(i);
    %insert constraints from blood conc + Ham's media
    modelGrowth = constrainMedium(inModel, bloodData.totMets, ux, false, true);
    modelGrowth.lb(strcmp(modelGrowth.rxns,'MAR03964')) = cellMaintenance;
    modelGrowth.lb(strcmp(modelGrowth.rxns,'f_MAR03964')) = cellMaintenance;
    modelGrowth.lb(strcmp(modelGrowth.rxns,'o_MAR03964')) = cellMaintenance;
    modelGrowth.c = double(strcmp(modelGrowth.rxns,'total_biomass'));%set objective function to cancer cell growth, including ECM
    %now, optimize for growth and see which metabolites that are exported
    %from fibroblasts to cancer cells
    %add those to the list and then block them by stopping the transport 
    %between f_s and s, and run again.
    %Then add the new metabolites to the list and so forth. Repeat until we
    %don't get any metabolites in the list.
    iteration = 1;
    
    %identify metabolites with infinite access
    %Filter mets that have >= 1000 as uptake bound
    %these are mets with reactions that have only one non-zero element in the S matrix
    %
    sComp = find(strcmp(modelGrowth.comps,'s'));
    sMetsSel = modelGrowth.metComps == sComp;
    openInExchRxns = ((sum(modelGrowth.S ~= 0, 1) == 1) & (sum(modelGrowth.S > 0, 1) == 1)).' & modelGrowth.ub >= 1000;
    openMetsSel = sMetsSel & sum(modelGrowth.S(:,openInExchRxns) ~= 0, 2) > 0;
    %TC001a
    %modelGrowth.metNames(openMetsSel) - looks good
    
    
    while true
        disp(['iteration: ' num2str(iteration)]);
        sol = solveLP(modelGrowth,1);
        if (sol.stat ~= 1)
           break; %succeeded will be false 
        end
        %1. Find the metabolites which are exported from the fibroblasts
        %and imported to the cancer cells
        %first find the reactions exporting these metabolites:
        %find the S compartment mets (use the input model as template)
        sComp = find(strcmp(modelGrowth.comps,'s'));
        sfComp = find(strcmp(modelGrowth.comps,'f_s'));
        soComp = find(strcmp(modelGrowth.comps,'o_s'));
        sMetsSel = modelGrowth.metComps == sComp;
        sfMetsSel = modelGrowth.metComps == sfComp;
        soMetsSel = modelGrowth.metComps == soComp;
        %now, find reactions with export to s from s_f
        rxnsExpFSelTmp = (sum(modelGrowth.S(sMetsSel,:) > 0,1) > 0).';
        rxnsExpFSelTmp2 = (sum(modelGrowth.S(sfMetsSel,:) < 0,1) > 0).'; %these are not just transport reactions
        %sum(rxnsExpFSelTmp)%4351
        %sum(rxnsExpFSelTmp2)%4154
        rxnsExpFSel = rxnsExpFSelTmp & rxnsExpFSelTmp2; %these are all export reactions from fibroblasts to S, 1064 in total
        %TC001b
        %constructEquations(modelGrowth, modelGrowth.rxns(rxnsExpFSel)) %looks good
        
        %filter out the reactions that don't carry flux
        actRxnsExpFSel = rxnsExpFSel & (sol.x ~= 0);%86 for example
        %Get the metabolites exported from the fibroblasts
        %Since all fluxes are positive and the transport between the s compartments is trivial,
        %we can simply sum up the fluxes of each metabolite over all reactions, and see if it is positive
        expMetsFSel = sMetsSel & (sum(modelGrowth.S(:,actRxnsExpFSel),2) > 0);
        %TC001c
        %sum(expMetsFSel) %14
        %modelGrowth.metNames(expMetsFSel)
        %constructEquations(modelGrowth, modelGrowth.rxns(actRxnsExpFSel))
        %looks good
        
        %Now, check which of those metabolites are taken up by the cancer
        %cells
        %A way to figure this out is to:
        %1. From the list of all reactions, remove all reactions that has a metabolite that
        %   belongs to the f_s or o_s compartments. In addition, remove those that only has 
        %   one metabolite. The idea is that after this filtering, no reactions will remain that
        %   moves any s metabolites anywhere except into the cancer cells.
        %2. We could then just sum up all fluxes for each s metabolite to see
        %   if it is taken up by the cancer cells. A negative sum would then indicate this.
        %This is easier than looking directly at import into the cancer cells, since those reactions are 
        %sometimes complicated, handling several metabolites, importing and exporting at the same time.
        boolS = modelGrowth.S ~= 0;
        rxnsNoFs = sum(boolS(sfMetsSel,:), 1) == 0; %all rxns except the ones using f_s mets
        rxnsNoOs = sum(boolS(soMetsSel,:), 1) == 0; %all rxns except the ones using o_s mets
        rxnsNoExch = sum(boolS, 1) ~= 1; %all rxns not having only one metabolite
        rxnSel = rxnsNoFs & rxnsNoOs & rxnsNoExch;
        %TC001d
        %check the reactions that use glucose[s]
        %glucMetInd = find(strcmp(modelGrowth.metNames, 'glucose') & sMetsSel);%2879
        %rxnsWithGluc = modelGrowth.S(glucMetInd,:) ~= 0;
        %constructEquations(modelGrowth,modelGrowth.rxns(rxnsWithGluc))%includes both glucose import and export + from f_s,o_s to s
        %constructEquations(modelGrowth,modelGrowth.rxns(rxnSel & rxnsWithGluc))%now these are gone, ok
        
        
        %sum(rxnSel) %34012, seems somewhat reasonable
        %now, get the metabolite fluxes. This is done as S*v, where the reactions are filtered
        metFluxes = modelGrowth.S(:,rxnSel) * sol.x(rxnSel);
        %now filter the metabolites to only include those in the s compartment
        metFluxesFilt = metFluxes;
        metFluxesFilt(~sMetsSel) = 0;
        impMetsSel = metFluxesFilt < 0;
        %TC001e
        %modelGrowth.metNames(impMetsSel)
        %These look reasonable, it is the ones from the diffusion model
        %sum(impMetsSel)%69
        %modelGrowth.metNames(expMetsFSel & impMetsSel)
        %modelGrowth.metNames(expMetsFSel & impMetsSel& ~openMetsSel) %nonzero at least
        
        %now, look at the mets that are exported from s and imported into c, skipping mets with unlimited access
        helpingMetsSel = expMetsFSel & impMetsSel & ~openMetsSel;
        modelGrowth.metNames(helpingMetsSel)
        
        %find the in exchange reactions for these mets, and check that the
        %input is larger than the upper bound of the exchange reaction
        %it seems this is not that important, only two metabolites, skip it
        %influxC = modelGrowth.S(:,actRxnsImpCSel) * sol.x(actRxnsImpCSel); %matrix multiplication
        %the exchange rxns:
        %exchRxns = (sum(modelGrowth.S ~= 0,1) == 1) & (sum(modelGrowth.S(helpingMetsSel,:) > 0,1) == 1);
        disp(['Found collaboration mets: ' num2str(sum(helpingMetsSel))])
        if sum(helpingMetsSel) == 0
            res.succeeded(i) = true;
            break;%we have found them all
        end
        %Add the mets as collaboration mets
        res.collaborationMets(:,i) = res.collaborationMets(:,i) | helpingMetsSel;
        
        %now, change the output of these metabolites from the fibroblasts to the  to the S2 compartment -
        %we assume those metabolites exist there, they should
        %do it with a loop to be safe
        collMetInd = find(helpingMetsSel);
        collMetNames = modelGrowth.metNames(collMetInd);
        %s2MetsSel = modelGrowth.metComps == sfComp;
        modelGrowth2 = modelGrowth; %remove this later
        %constructEquations(modelGrowth, modelGrowth.rxns(actRxnsExpFSel))
        for j = 1:length(collMetInd)
            sfMetSel = strcmp(modelGrowth.metNames,collMetNames{j}) & sfMetsSel;
            if (sum(sfMetSel) ~= 1)
                disp('something went wrong')
            end
            sel2 = rxnsExpFSel & (modelGrowth2.S(collMetInd(j),:) > 0).'; 
            %constructEquations(modelGrowth2, modelGrowth2.rxns(sel2)) %{'(13Z)-eicosenoic acid[f_s] => (13Z)-eicosenoic acid[s]'}
            %modelGrowth2.S(s2MetSel,sel2) = modelGrowth2.S(collMetInd(j),sel2);
            modelGrowth2.S(collMetInd(j),sel2) = 0;
            %constructEquations(modelGrowth2, modelGrowth2.rxns(sel2)) %{'(13Z)-eicosenoic acid[f_s] => '}, ok
        end
        %constructEquations(modelGrowth2, modelGrowth2.rxns(actRxnsExpFSel)) %looks good
        modelGrowth = modelGrowth2;
        iteration = iteration + 1;
    end
end

