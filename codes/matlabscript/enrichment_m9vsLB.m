load('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/M9andLB_DErxns.mat')
load('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/iML1515.mat')
model=iML1515;
common_rxns=intersect(LBDEitem,M9DEitem);

out_file=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/','M9_LB_common_rxns_enrichment.txt')
DEflag=2;

DEitem=common_rxns;

subSystemEnrichment(model,DEitem,out_file,DEflag)


common_rxns=setdiff(LBDEitem,M9DEitem);

out_file=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/','only_LB_enrichment.txt')
DEflag=2;

DEitem=common_rxns;

subSystemEnrichment(model,DEitem,out_file,DEflag)



common_rxns=setdiff(M9DEitem,LBDEitem);

out_file=strcat('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/','only_M9_enrichment.txt')
DEflag=2;

DEitem=common_rxns;

subSystemEnrichment(model,DEitem,out_file,DEflag)

