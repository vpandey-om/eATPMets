%% plot flux for M9

% 
% rxns={'ACGAMT';'ACGK';'ACGS';'ADEt2rpp';'AHCYSNS';'AMPN';'CYTDK2';'DADA';'NNAM'; ...
%     'PAPPT3';'PUNP2';'PUNP6';'RPI';'TKT1';'XPPT'}
rxns={'ACGK';'ACGS';'AMPN';'CYTDK2';'DADA';'NNAM';'PUNP2';'PUNP6';'TKT1';'GLUPRT';'GTHS'}

rxns={'ACGK';'ACGS';'AMPN';'CYTDK2';'DADA';'PUNP2';'PUNP6';'TKT1';'GTHS'}

load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9_alternative_final.mat')
NFind=getAllVar(model,{'NF'});
PNFind=getAllVar(model,{'PERTURB_NF'});
ctl_fluxes=sol_matrix(NFind,1:20);
eatp_fluxes=sol_matrix(PNFind,1:20);
[~,ind]=ismember(rxns,model.rxns);

sol1Arr=ctl_fluxes(ind,:);
sol2Arr=eatp_fluxes(ind,:);

condi=[];
wtFlux=[];
compound=[];
for i=1:numel(rxns)

    for j=1:size(sol1Arr,2)
        condi=[condi;{'M9'}];
        wtFlux=[wtFlux;sol1Arr(i,j)];
        compound=[compound;{rxns{i}}];
    end


    for j=1:size(sol2Arr,2)
        condi=[condi;{'M9+ATP'}];
        wtFlux=[wtFlux;sol2Arr(i,j)];
        compound=[compound;{rxns{i}}];
    end
    
    mltable = table(compound,wtFlux,condi);

% boxchart(mltable.compound,mltable.wtFlux,'GroupByColor',mltable.condi)
% ylabel('Maximum production rate')
% legend

end
filename=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/m9_rxns_flux_plot.txt')
writetable(mltable,filename,'Delimiter','tab')


%% lb plot


%% plot flux for M9

% 
% rxns={'ACGAMT';'ACGK';'ACGS';'ADEt2rpp';'AHCYSNS';'AMPN';'CYTDK2';'DADA';'NNAM'; ...
%     'PAPPT3';'PUNP2';'PUNP6';'RPI';'TKT1';'XPPT'}
rxns={'GALUi','UMPK','INSH','TRPtipp'};
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/LB_alternative_final.mat')
NFind=getAllVar(model,{'NF'});
PNFind=getAllVar(model,{'PERTURB_NF'});
ctl_fluxes=sol_matrix(NFind,1:20);
eatp_fluxes=sol_matrix(PNFind,1:20);
[~,ind]=ismember(rxns,model.rxns);

sol1Arr=ctl_fluxes(ind,:);
sol2Arr=eatp_fluxes(ind,:);

condi=[];
wtFlux=[];
compound=[];
for i=1:numel(rxns)

    for j=1:size(sol1Arr,2)
        condi=[condi;{'LB'}];
        wtFlux=[wtFlux;sol1Arr(i,j)];
        compound=[compound;{rxns{i}}];
    end


    for j=1:size(sol2Arr,2)
        condi=[condi;{'LB+ATP'}];
        wtFlux=[wtFlux;sol2Arr(i,j)];
        compound=[compound;{rxns{i}}];
    end
    
    mltable = table(compound,wtFlux,condi);

% boxchart(mltable.compound,mltable.wtFlux,'GroupByColor',mltable.condi)
% ylabel('Maximum production rate')
% legend

end
filename=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/lb_rxns_flux_plot.txt')
writetable(mltable,filename,'Delimiter','tab')








% boxchart(mltable.compound,mltable.wtFlux,'GroupByColor',mltable.condi)
% ylabel('Maximum production rate')
% legend

    subsys=agaDraftModel.subSystems;
    for i=1:numel(subsys)
        if(numel(agaDraftModel.subSystems{i})>1)
        subsys{i}=strjoin(agaDraftModel.subSystems{i},{'|'});
        end

    end
