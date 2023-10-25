load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9_dereg_mets.mat')
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/iML1515.mat')
% read tables
m9=readtable('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/M9_rxn_pvals_kstest.txt','Delimiter','tab');
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9datamodels.mat')

cofactors={'CO2 CO2','Copper','Calcium','H+','Nitrite','Fluoride','Zinc', 'Co2+','Ammonium','Potassium','Cu+','Fe2+ mitochondria','H2O H2O', ...
'ATP C10H12N5O13P3','Phosphate','ADP C10H12N5O10P2','Nicotinamide adenine dinucleotide','Diphosphate','Nicotinamide adenine dinucleotide - reduced', ...
'Nicotinamide adenine dinucleotide phosphate','Nicotinamide adenine dinucleotide phosphate - reduced','Coenzyme A', ...
'AMP C10H12N5O7P','Acyl carrier protein','Phosphate','O2 O2','CMP C9H12N3O8P','GTP C10H12N5O14P3','CTP C9H12N3O14P3', ...
'UDP C9H11N2O12P2','Sodium','GDP C10H12N5O11P2','FMN C17H19N4O9P','Sulfite'};
cofator_mets={};
for i=1:numel(cofactors)
    ind=find(ismember(iML1515.metNames,cofactors{i}));
    cofator_mets=[cofator_mets;iML1515.mets(ind)]
end

cofator_mets=unique(cofator_mets);

metData=M9deregmodelmets;
mets=cellstr(metData(2:end,1));
metNames=cellstr(metData(2:end,2));
FC=metData(2:end,5);
FC=cell2mat(FC);
cytoMets={};
for i=1:numel(mets)
    if strcmp('c',mets{i}(end))
        cytoMets=[cytoMets;mets{i}];
    end
end

cytoMets=strrep(cytoMets,'M_','');

totalDEGenes=model_genes(model_gene_ratios<2/3 | model_gene_ratios>3/2);

result1={};
result2={};
result3={};
for i=1:numel(cytoMets)
    
    [~,ind]=ismember(cytoMets{i},iML1515.mets);    
    rxns=iML1515.rxns(find(iML1515.S(ind,:)));
    derxns=intersect(m9.DEitem,rxns);
    deFC=intersect(rxns,dereg_rxns);
    [~,ind]=ismember(rxns,iML1515.rxns);
    [rr,gg,~]=find(iML1515.rxnGeneMat(ind,:));
    
    % for de genes 
    
    deGenes=intersect(iML1515.genes(unique(gg)),totalDEGenes);
    result1{i}=[derxns;deGenes];
    
    %% round2
    
    % find metabolites 
    rxns=iML1515.rxns(unique(rr));
     [~,ind]=ismember(rxns,iML1515.rxns);
    [rr,mm,~]=find(iML1515.S(:,ind));
    
    m=setdiff(iML1515.mets(unique(mm)),cofator_mets);
    [~,ind]=ismember(m,iML1515.mets);
    [rr,mm,~]=find(iML1515.S(ind,:));
    
    rxns=iML1515.rxns(unique(rr));
    derxns=intersect(m9.DEitem,rxns);
    deFC=intersect(rxns,dereg_rxns);
    [rr,gg,~]=find(iML1515.rxnGeneMat(unique(rr),:));
    deGenes=intersect(iML1515.genes(unique(gg)),totalDEGenes);
    result2{i}=[derxns;deGenes];
   
    
    
    
    
    
    
end

%%

allitems={};
for i=1:numel(result1)
    if ~(i==11)
        allitems=[allitems;result1{i};result2{i};cytoMets{i}];  
    end
    
end
allitems=unique(allitems)


%% For LB data 
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/LB_deReg_mets.mat')
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/iML1515.mat')
% read tables
lb=readtable('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/LB_rxn_pvals_kstest.txt','Delimiter','tab');
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/LBdatamodels.mat')

cofactors={'CO2 CO2','Copper','Calcium','H+','Nitrite','Fluoride','Zinc', 'Co2+','Ammonium','Potassium','Cu+','Fe2+ mitochondria','H2O H2O', ...
'ATP C10H12N5O13P3','Phosphate','ADP C10H12N5O10P2','Nicotinamide adenine dinucleotide','Diphosphate','Nicotinamide adenine dinucleotide - reduced', ...
'Nicotinamide adenine dinucleotide phosphate','Nicotinamide adenine dinucleotide phosphate - reduced','Coenzyme A', ...
'AMP C10H12N5O7P','Acyl carrier protein','Phosphate','O2 O2','CMP C9H12N3O8P','GTP C10H12N5O14P3','CTP C9H12N3O14P3', ...
'UDP C9H11N2O12P2','Sodium','GDP C10H12N5O11P2','FMN C17H19N4O9P','Sulfite'};
cofator_mets={};
for i=1:numel(cofactors)
    ind=find(ismember(iML1515.metNames,cofactors{i}));
    cofator_mets=[cofator_mets;iML1515.mets(ind)]
end

cofator_mets=unique(cofator_mets);

metData=richderegmodelmets2;
mets=cellstr(metData(2:end,1));
metNames=cellstr(metData(2:end,2));
FC=metData(2:end,5);
FC=cell2mat(FC);
cytoMets={};
for i=1:numel(mets)
    if strcmp('c',mets{i}(end))
        cytoMets=[cytoMets;mets{i}];
    end
end

cytoMets=strrep(cytoMets,'M_','');

totalDEGenes=model_genes(model_gene_ratios<2/3 | model_gene_ratios>3/2);

result1={};
result2={};
result3={};
for i=1:numel(cytoMets)
    
    [~,ind]=ismember(cytoMets{i},iML1515.mets);    
    rxns=iML1515.rxns(find(iML1515.S(ind,:)));
    derxns=intersect(lb.DEitem,rxns);
    deFC=intersect(rxns,dereg_rxns);
    [~,ind]=ismember(rxns,iML1515.rxns);
    [rr,gg,~]=find(iML1515.rxnGeneMat(ind,:));
    
    % for de genes 
    
    deGenes=intersect(iML1515.genes(unique(gg)),totalDEGenes);
    result1{i}=[derxns;deGenes];
    
    % round2
    
    % find metabolites 
    rxns=iML1515.rxns(unique(rr));
     [~,ind]=ismember(rxns,iML1515.rxns);
    [rr,mm,~]=find(iML1515.S(:,ind));
    
    m=setdiff(iML1515.mets(unique(mm)),cofator_mets);
    [~,ind]=ismember(m,iML1515.mets);
    [rr,mm,~]=find(iML1515.S(ind,:));
    
    rxns=iML1515.rxns(unique(rr));
    derxns=intersect(lb.DEitem,rxns);
    deFC=intersect(rxns,dereg_rxns);
    [rr,gg,~]=find(iML1515.rxnGeneMat(unique(rr),:));
    deGenes=intersect(iML1515.genes(unique(gg)),totalDEGenes);
    result2{i}=[derxns;deGenes];
   
    
    
    
    
    
    
end

%%

result1{5}={'INSH';'PUNP1';'URIH'};

allitems={};
for i=1:numel(result1)
    if ~(i==11)
        allitems=[allitems;result1{i};result2{i};cytoMets{i}];  
    end
    
end
allitems=unique(allitems)

% [~,ind]=ismember({'b0030';'b4384'},iML1515.genes); {'INSH','PUNP1'}