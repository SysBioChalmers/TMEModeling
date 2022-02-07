%Adds a small protein usage cost (on the cancer cells) to the transport reaction from f_s to s,
%to avoid random collaborations. This will lead to that the cell will prefer 
%not to use collaboration mets if it doesn't have to to maximize growth.
function outModel = addCollaborationCost(inModel)
outModel = inModel;
sComp = find(strcmp(outModel.comps,'s'));
sfComp = find(strcmp(outModel.comps,'f_s'));
sMetsSel = outModel.metComps == sComp;
sfMetsSel = outModel.metComps == sfComp;
rxnsExpFSelTmp = (sum(outModel.S(sMetsSel,:) > 0,1) > 0).';
rxnsExpFSelTmp2 = (sum(outModel.S(sfMetsSel,:) < 0,1) > 0).'; 
rxnsExpFSel = rxnsExpFSelTmp & rxnsExpFSelTmp2; %these are all export reactions from fibroblasts to S, 1064 in total
%constructEquations(outModel, outModel.rxns(rxnsExpFSel))
protPoolSel = strcmp(outModel.mets, 'prot_pool');
outModel.S(protPoolSel,rxnsExpFSel) = -10^-6;
%constructEquations(outModel, outModel.rxns(rxnsExpFSel)) %looks good

end

