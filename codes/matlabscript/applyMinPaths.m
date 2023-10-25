params.timeLimit=1800;
params.cores=4;
model=iML1515;
model=removeRxns(model,{'BIOMASS_Ec_iML1515_core_75p37M','BIOMASS_Ec_iML1515_WT_75p37M'});
[finalVectors] = detMinSpan(model, params);