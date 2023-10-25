% growth analysis

load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9_alternative_final.mat')
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9modelwithData.mat')
tmpModel=coMExpMet;
% fit growth
[~,bind]=ismember('F_BIOMASS_Ec_iML1515_core_75p37M',tmpModel.varNames);
[~,pbind]=ismember('PERTURB_F_BIOMASS_Ec_iML1515_core_75p37M',tmpModel.varNames);

tmpModel.var_lb(bind)=0;
tmpModel.var_lb(pbind)=0;
sol1Arr=[];
sol2Arr=[];

for i=20
 tmp=tmpModel;
 soln=sol_matrix(:,i);
 binVar=[tmpModel.relExp.forB;tmpModel.metB];
 alwayson=binVar(soln(binVar)>0.98);
 alwaysoff=binVar(soln(binVar)<0.08);
 tmp.var_lb(alwayson)=1;
 tmp.var_ub(alwaysoff)=0;

%  saveRes=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/bbb/BBB_M9_alt', num2str(i),'.mat')
%     
%  [BBB, BBBprod1, BBBprod2] = testBBBinREMI_time(tmp,300,saveRes);
    
 
end

%% load BBBs 
num=20;
BBB1=[];
BBB2=[];
for i=1:num
    
    load(strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/bbb/BBB_M9_alt', num2str(i),'.mat'))
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
BBBr_name{2}='Riboflavin';
BBBr_name{5}='FAD oxidized';

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
writetable(mltable,'/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/M9_for_box_plot.txt','Delimiter','tab')

result1.ratios=ratios;
result1.BBBr_name=BBBr_name;
result1.BBBr=BBBr;
result1.BBB=BBB;
