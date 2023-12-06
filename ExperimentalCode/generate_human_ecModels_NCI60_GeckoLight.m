% This script requires you to download the zenodo from "An atlas of human metabolism", https://zenodo.org/record/3583004 
% This is a modified version of generate_human_ecModels_NCI60, where we generate GECKO models using GeckoLight instead of the original GECKO.
%

%We write all files in the Human-GEM zenodo file structure
%for practical reasons, we add our scripts to the paths, there are a lot of cd etc. below
addpath('C:/Work/MatlabCode/projects/TMEModeling/TMEModeling/GeckoLightEval') %Don't save the path, it is only temporarily added
cd C:/Work/MatlabCode/projects/HumanGemZenodo/Human1_Publication_Data_Scripts/ec_GEMs/ComplementaryScripts

clear %this is important since we use who below
load('../models/humanGEM_cellLines/11models.mat')
modelNames = who;
current    = pwd;
%Generate enzyme-constrained models
for i=1:length(modelNames)
    disp(i)
    cd (current)
    cellName = modelNames{i};
    mkdir (['../models/' cellName '_gl'])
    cd (['../models/' cellName  '_gl'])
    save([cellName '_gl.mat'],cellName)
    mkdir Data
    cd (current)
    %Models mat files are saved in their respective folder by enhanceGEM_cellLine
    enhanceGEM_cellLine_GeckoLight(cellName);
end
cd (current)
