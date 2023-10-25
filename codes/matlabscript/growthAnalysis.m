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

for i=1:20
 tmp=tmpModel;
 soln=sol_matrix(:,i);
 binVar=[tmpModel.relExp.forB;tmpModel.metB];
 alwayson=binVar(soln(binVar)>0.98);
 alwaysoff=binVar(soln(binVar)<0.08);
 tmp.var_lb(alwayson)=1;
 tmp.var_ub(alwaysoff)=0;
 %%% 
 f=zeros(numel(tmp.f),1);
 f(bind)=1;
 tmp.f=f;
 sol=solveTFAmodel(tmp);
 sol1Arr(i)=sol.val;
 
 %%% 
 f=zeros(numel(tmp.f),1);
 f(pbind)=1;
 tmp.f=f;
 sol=solveTFAmodel(tmp);
 sol2Arr(i)=sol.val;
 
end



save('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9_growth_rate', 'sol1Arr','sol2Arr');
%% growth analysis for LB data 

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
 %%% 
 f=zeros(numel(tmp.f),1);
 f(bind)=1;
 tmp.f=f;
 sol=solveTFAmodel(tmp);
 sol1Arr(i)=sol.val;
 
 %%% 
 f=zeros(numel(tmp.f),1);
 f(pbind)=1;
 tmp.f=f;
 sol=solveTFAmodel(tmp);
 sol2Arr(i)=sol.val;
 
end

 save('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/LB_growth_rate', 'sol1Arr','sol2Arr');




