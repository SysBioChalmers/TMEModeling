%Changes for example glucose[f_e] => glucose[e] to glucose[f_e] => 
function outModel = blockCollaboration(inModel, metsToKeep)
if nargin < 2
    metsToKeep = {};
end
outModel = inModel;
        %1. Find the metabolites which are exported from the fibroblasts
        %and imported to the cancer cells
        %first find the reactions exporting these metabolites:
        %find the S compartment mets (use the input model as template)
sComp = find(strcmp(outModel.comps,'e'));
sfComp = find(strcmp(outModel.comps,'f_e'));
sMetsSel = outModel.metComps == sComp;
sfMetsSel = outModel.metComps == sfComp;
%Filter out the mets to keep
metsToKeepSel = ismember(outModel.metNames, metsToKeep);

rxnsExpFSelTmp = (sum(outModel.S(sMetsSel & ~metsToKeepSel,:) > 0,1) > 0).';
rxnsExpFSelTmp2 = (sum(outModel.S(sfMetsSel  & ~metsToKeepSel,:) < 0,1) > 0).'; 
rxnsExpFSel = rxnsExpFSelTmp & rxnsExpFSelTmp2; %these are all export reactions from fibroblasts to S, 1064 in total


outModel.S(sMetsSel,rxnsExpFSel) = 0;

end

