% growth analysis

load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/LB_alternative_final.mat')
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/LBmodelwithData.mat')
tmpModel=coMExpMet;
% fit growth
[~,bind]=ismember('F_BIOMASS_Ec_iML1515_core_75p37M',tmpModel.varNames);
[~,pbind]=ismember('PERTURB_F_BIOMASS_Ec_iML1515_core_75p37M',tmpModel.varNames);

tmpModel.var_lb(bind)=0;
tmpModel.var_lb(pbind)=0;
sol1Arr=[];
sol2Arr=[];

for i=1:20
 tmp=tmpModel;
 soln=sol_matrix(:,i);
 binVar=[tmpModel.relExp.forB;tmpModel.metB];
 alwayson=binVar(soln(binVar)>0.98);
 alwaysoff=binVar(soln(binVar)<0.08);
 tmp.var_lb(alwayson)=1;
 tmp.var_ub(alwaysoff)=0;

%  saveRes=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/bbb/BBB_LB_alt', num2str(i),'.mat')
%     
%  [BBB, BBBprod1, BBBprod2] = testBBBinREMI_time(tmp,300,saveRes);
    
 
end

%% load BBBs 
num=20;
BBB1=[];
BBB2=[];
for i=1:num
    
    load(strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/bbb/BBB_LB_alt', num2str(i),'.mat'))
    BBB1=[BBB1 BBBprod1];
    BBB2=[BBB2 BBBprod2];
    
end

BBB1=cell2mat(BBB1);
BBB2=cell2mat(BBB2);

% calculate geometric mean
ratios=BBB1;
for i=1:size(BBB1,1)
    for j=1:size(BBB1,2)
        ratios(i,j)=BBB2(i,j)/BBB1(i,j);
    end
end
allBBB=[BBB1 BBB2];
georatios = geomean(ratios,2);

%%
deBBB=allBBB(georatios>1.5 | georatios<(2/3),:);


condi={};
BBBr=BBB(georatios>1.5 | georatios<(2/3) );
[~,ind]=ismember(BBBr,coMExpMet.mets);
BBBr_name=coMExpMet.metNames(ind);


wtFlux=[];
compound=[];

for i=1:size(deBBB,1)
    for j=1:num
        condi=[condi;{'control'}];
        wtFlux=[wtFlux;deBBB(i,j)];
        compound=[compound;BBBr_name(i)];       
    end
    for j=num+1:2*num
        condi=[condi;{'eATP'}];
        wtFlux=[wtFlux;deBBB(i,j)];
        compound=[compound;BBBr_name(i)];
    end
end

mltable = table(compound,wtFlux,condi);

% boxchart(mltable.compound,mltable.wtFlux,'GroupByColor',mltable.condi)
% ylabel('Maximum production rate')
% legend
writetable(mltable,'/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/LB_for_box_plot.txt','Delimiter','tab')

result2.ratios=ratios;
result2.BBBr_name=BBBr_name;
result2.BBBr=BBBr;
result2.BBB=BBB;

%% combine result1 and result 2
BBBc=union(result1.BBBr,result2.BBBr);
[~,ind]=ismember(BBBc,coMExpMet.mets);
BBBc_name=coMExpMet.metNames(ind);
BBBc_name{7}='Riboflavin';
BBBc_name{2}='FAD oxidized';
[BBBc BBBc_name]
condi=[]
wtFlux=[];
compound=[];
[~,ind]=ismember(BBBc,result1.BBB);
deBBB=[result1.ratios(ind,:) result2.ratios(ind,:)];
for i=1:size(deBBB,1)
    for j=1:num
        condi=[condi;{'M9'}];
        wtFlux=[wtFlux;log2(deBBB(i,j))];
        compound=[compound;BBBc_name(i)];       
    end
    for j=num+1:2*num
        condi=[condi;{'LB'}];
        wtFlux=[wtFlux;log2(deBBB(i,j))];
        compound=[compound;BBBc_name(i)];
    end
end

mltable = table(compound,wtFlux,condi);

% boxchart(mltable.compound,mltable.wtFlux,'GroupByColor',mltable.condi)
% ylabel('Maximum production rate')
% legend
writetable(mltable,'/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/LB_M9_for_box_plot.txt','Delimiter','tab')



