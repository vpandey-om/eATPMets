load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9_alternative_final.mat')
NFind=getAllVar(model,{'NF'});
PNFind=getAllVar(model,{'PERTURB_NF'});
ctl_fluxes=sol_matrix(NFind,1:20);
eatp_fluxes=sol_matrix(PNFind,1:20);

ctl_mean_flux=mean(ctl_fluxes,2);
eatp_mean_flux=mean(eatp_fluxes,2);
ctl_var_flux=var(ctl_fluxes,0,2);
eatp_var_flux=var(eatp_fluxes,0,2);
%% calculate z score for diffrential analysis 
zscore_flux=(eatp_mean_flux-ctl_mean_flux)./sqrt(ctl_var_flux+eatp_var_flux);

pvals=ones(size(ctl_fluxes,1),1);
for i=1:size(ctl_fluxes,1)
    p = normcdf(zscore_flux(i,:));
%     p = signrank(ctl_fluxes(i,:),eatp_fluxes(i,:));
    pvals(i)=p;
end

%%
%%% get based on z score 

cut=1;
see=[ctl_fluxes(abs(zscore_flux)>cut,:) eatp_fluxes(abs(zscore_flux)>cut,:)];

xx=find(abs(zscore_flux)>cut);

boxplot([ctl_fluxes(xx(3),:)',eatp_fluxes(xx(3),:)'])

pvals=[]
for i=1:numel(xx)
    pvals(i)= signrank(ctl_fluxes(xx(i),:),eatp_fluxes(xx(i),:));
end


%%
[uE,~,ix] = unique(model.subSystems(xx));
tally = accumarray(ix, 1)/length(xx);
Out = table(uE, tally, 'VariableNames',{'Subsystems','Frequency'})

%% rxn enrichment analysis 
out_file=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/','M9_rxns_enrichment.txt')
DEflag=2;

DEitem=model.rxns(xx);

subSystemEnrichment(model,DEitem,out_file,DEflag)

writematrix([ctl_fluxes eatp_fluxes],'/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/M9_flux_data.txt','Delimiter','tab')
pvals=pvals';

% T = table(DEitem,pvals)
subsys=model.subSystems(xx);
score=zscore_flux(xx);

T = table(DEitem,subsys,pvals,score)

writetable(T,'/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/M9_rxn_pvals_kstest.txt','Delimiter','tab');


% pvals=ones(size(ctrl_flux,1),1);
% for i=1:size(ctrl_flux,1)
% [h,p] = kstest2(ctrl_flux(i,:),perb_flux(i,:),'Alpha',0.05);
% pvals(i)=p;
% end
%%
txx=[ctl_fluxes eatp_fluxes];
xx=txx;
xx(xx<-1e-8 | xx >1e-8)=1;
inactiveRxns=model.rxns(sum(xx,2)==0);



