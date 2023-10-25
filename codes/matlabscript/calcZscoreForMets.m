function [z_score,mets,store_size]=calcZscoreForMets(model,zRxns)
     %% model is metabolic model in cobra format
     % pvals is pvalues for the each enzyme based on t test     
     %cofac_file is a excel file of cofactors. For these cofactors
     %analysis will not take place.
     mets=model.mets;
     store_size=NaN(numel(mets),1);
     z_score=NaN(numel(mets),1);
     for i=1:numel(mets)
         [crap,ind]=ismember(mets{i},model.mets);
         rxnInd=find(model.S(ind,:));
         met_score=zRxns(rxnInd);
         met_score=met_score(~isnan(met_score));
         if numel(met_score)>0
            store_size(i)=numel(met_score);
            z_score(i)=sum(met_score)*(1/sqrt(numel(met_score)));
         end
         
     end
      
end
    