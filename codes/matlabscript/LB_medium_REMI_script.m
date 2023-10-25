% get parent folder from current folder 
eatpMetabolism = cd;
modelFolder=strcat(eatpMetabolism,'/models');
scriptFolder=strcat(eatpMetabolism,'/scripts');
dataFolder=strcat(eatpMetabolism,'/data');
matFolder=strcat(eatpMetabolism,'/matfiles');
addpath(modelFolder)
addpath(scriptFolder)

%% ask folder for setting path for software 
%remipath = input('please provide REMI folder path: (/Users/vikash/Documents/MATLAB/matpacks/REMI2.0):',"s")
remipath='/Users/vikash/Documents/MATLAB/matpacks/REMI2.0'
addpath(genpath(remipath)) % REMI method 

%matTFApath = input('please provide matTFA folder path: (/Users/vikash/Documents/MATLAB/matpacks/matTFA):',"s")
matTFApath='/Users/vikash/Documents/MATLAB/matpacks/matTFA'
addpath(genpath(matTFApath)) % matTFA

%texFBApath = input('please provide texFBA folder path: (/Users/vikash/Documents/MATLAB/matsoftware/FBA_Toolboxes/TranscriptomicsIntegration):',"s")
texFBApath='/Users/vikash/Documents/MATLAB/matsoftware/FBA_Toolboxes/TranscriptomicsIntegration'
addpath(genpath(texFBApath)) % TranscriptomicsIntegration



%% load models and metabolite data 
% apply REMI on one metabolite data sets 
load(fullfile(modelFolder,'LBiML1515.mat'))
iML1515=model;
% get match experimental metabolite data to models 
load(fullfile(matFolder,'LB_deReg_mets.mat'))
metData=richderegmodelmets2;
mets=cellstr(metData(2:end,1));
FC=metData(2:end,5);
FC=cell2mat(FC);
% compsymbol=cellfun(@(x) x(end), mets, 'UniformOutput', false);

% remove extracellular metabolites 
intraMets={};
itraFC=[];
for i=1:numel(mets)
    if ~ strcmp('e',mets{i}(end))
        intraMets=[intraMets;mets{i}];
        itraFC=[itraFC;FC(i)];
    end
    
end


dereg_metFC=itraFC;
dereg_mets=intraMets;


%% get expression data deregualted one.
load(fullfile(matFolder,'LB_genes.mat'))
geneData=Allgfpvfinal04ST2018reduced;
% matches = strfind(cellstr(geneData(:,7)),',');
% tf = any(vertcat(matches{:}))
% find match genes with model
% geneExp=cellstr(geneData(2:1699,7));
% geneLC2fold=cell2mat(geneData(2:1699,39));
geneExp=cellstr(geneData(2:end,7));
tmpdata=geneData(2:end,39);
idx=[]
for i=1:numel(tmpdata)
    if isa(tmpdata{i},'double')
        idx=[idx;i];
        
    end
end
geneExp=geneExp(idx);
geneLC2fold=cell2mat(tmpdata(idx));
geneRatio=power(2,geneLC2fold);
[~,ind]=ismember(iML1515.genes,geneExp);
model_gene_data=[iML1515.genes(ind>0) geneExp(ind(ind>0)) num2cell(geneRatio(ind(ind>0)))]
model_genes=iML1515.genes(ind>0);
model_gene_ratios=geneRatio(ind(ind>0));
%% apply REMI to get reaction expression
[rxns,rxnExp]=evaluateGPR(iML1515,model_genes,model_gene_ratios,@geomean,@mean);
% find is nan reaction expression
boolVal=~isnan(rxnExp);
rxns_all=rxns(boolVal)';
rxns_FC=rxnExp(boolVal);
% deregulated reactions
cut=1.5;
boolcond=(rxns_FC>cut) | (rxns_FC<(1/cut));
dereg_rxnFC=rxns_FC(boolcond);
dereg_rxns=rxns_all(boolcond);
% remove regualtion greater than 100 fold
dereg_rxnFC(dereg_rxnFC>100)=100;

%% make model in TFA format
iML1515.rev=ones(numel(iML1515.rev),1);
model=addUseVariablesDH(iML1515);
model=addNetFluxVariables(model);
sol=solveTFAmodel(model);
% put 90% of biomass as constraint
[~,biomassIdx]=ismember('F_BIOMASS_Ec_iML1515_core_75p37M',model.varNames);
model.var_lb(biomassIdx)=sol.val*0.9;

%% add expression constriant
netM=model;
coMexp=addRelConsExpression(netM,netM,dereg_rxns,dereg_rxnFC)
sol=solveTFAmodel(coMexp);

%% add relative metaboliets 
regMets=strrep(dereg_mets,'M_','');
regMetRatio=dereg_metFC;
[coMExpMet]=addRelMetabolite(coMexp,coMexp,regMets,regMetRatio,true);
sol2=solveTFAmodel(coMExpMet);
save(fullfile(matFolder,'LBdatamodels'),'regMets','regMetRatio','model_genes','model_gene_ratios','dereg_rxns','dereg_rxnFC')

%% alternative solutions
comIdx=[coMExpMet.relExp.forB;coMExpMet.metB];
path_save=fullfile(matFolder,'LB_alternative.mat')

addpath(genpath('/Users/vikash/Documents/MATLAB/matsoftware/FBA_Toolboxes/tFBA_vDec2014'))
findAltCombi(499,coMExpMet,comIdx,path_save,600);





