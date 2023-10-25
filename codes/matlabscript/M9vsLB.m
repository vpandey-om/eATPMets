%% Compare M9 and LB medium 


%%
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/iML1515.mat')
model=iML1515;
[secbool,upbool]=findExcRxns(model);
[~,~,LB_data]=xlsread('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/LB_medium.xlsx');
addpath(genpath('/Users/vikash/Documents/MATLAB/matpacks/REMI2.0')) % REMI method 
%%
LB_mets=LB_data(2:55,22);
LB_met_vals=str2double(LB_data(2:55,23));
LB_mets=strrep(LB_mets,')','');
LB_mets=strrep(LB_mets,'(','_');
[~,ind]=ismember(LB_mets,model.rxns);
[model.lb(ind) model.ub(ind)]

model.lb(ind)=-LB_met_vals;

 sol1=solveFBACplex(model)

 save('LBiML1515.mat', 'model')

%% M9 medium 
[~,~,M9_data]=xlsread('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/M9_medium.xlsx','M9');
M9_mets=M9_data(:,1);
M9_met_vals=str2double(M9_data(:,3));
model=iML1515;
sol=solveFBACplex(model)
% fix glucose
[~,ind]=ismember('EX_glc__D_e',model.rxns);

[~,ind]=ismember(M9_mets,model.rxns);
vals=[];
model.lb(ind)=M9_met_vals
sol2=solveFBACplex(model);
vals=[vals;sol2.objval];  
save('M9iML1515.mat', 'model')


