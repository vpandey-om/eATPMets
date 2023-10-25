% reporter metabolite analysis 
%% set path
eatpMetabolism = cd;
modelFolder=strcat(eatpMetabolism,'/models');
dataFolder=strcat(eatpMetabolism,'/data');
scriptFolder=strcat(eatpMetabolism,'/MatlabScripts');
matFolder=strcat(eatpMetabolism,'/matfiles');
addpath('/Users/vikash/Documents/MATLAB/matsoftware/Reporter_script') % for reporter script 
addpath(genpath(modelFolder))
addpath(genpath(scriptFolder))
addpath(genpath(matFolder))
addpath(genpath(dataFolder))

%% load  models 
load('iML1515.mat')
cofators={'CO2 CO2','O2 O2','H2O H2O','Nicotinamide adenine dinucleotide','ATP C10H12N5O13P3',...
    'H+','ADP C10H12N5O10P2','Phosphate','Coenzyme A','Diphosphate',...
    'GDP C10H12N5O11P2', 'GTP C10H12N5O14P3','UTP C9H11N2O15P3','Nicotinamide adenine dinucleotide - reduced',...
    'Nicotinamide adenine dinucleotide phosphate','Nicotinamide adenine dinucleotide phosphate - reduced'}
    
mets=[]
for i=1:numel(cofators)
  ind=find(ismember(iML1515.metNames,cofators{i}));
  mets=[mets;iML1515.mets(ind)];
end

%readdata 
load('lb_rxns_fold_change.mat')
rxns=lb_rxn_fold(:,1);
rxns_log2Fold=cell2mat(lb_rxn_fold(:,2));
rxns_zscore=(rxns_log2Fold-mean(rxns_log2Fold))./std(rxns_log2Fold);

%% apply reporter metaboliets 
rxnScore=nan(numel(iML1515.rxns),1);
[~,ind]=ismember(rxns,iML1515.rxns);
rxnScore(ind)=rxns_log2Fold;
model=iML1515;
% remove cofatores
[~,metIdx]=ismember(mets,model.mets);
model.S(metIdx,:)=0;

%% apply wit cobra 
nRand=10000;
pValFlag=false,
nLayers=1;
[normScore,nRxnsMet,nRxnsMetUni,rawScore] = reporterMets(model,rxnScore,nRand,pValFlag,nLayers);
%%
[z_score_mets,metids,store_size]=calcZscoreForMets(model,rxnScore);
cluster_num=unique(store_size(~isnan(store_size)));
zscore_without_nan=z_score_mets(~isnan(z_score_mets));
[muk,sigmak,matrix_zscore]=calRandZscoreForMets(zscore_without_nan,cluster_num,10000);

[z_correct]=correctedZscore(z_score_mets,store_size,cluster_num,muk,sigmak);

zCorrect_without_nan=z_correct(~isnan(z_correct));
mets=model.mets(~isnan(z_correct));
percentValue=prctile(z_correct,[80]);
ind=find(z_correct>percentValue);
%% calculate pvalues
[reporter_pvals]=estimatePvalues(z_score_mets,store_size,cluster_num,matrix_zscore);



