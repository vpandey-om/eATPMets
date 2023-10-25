function [muk,sigmak,matrix_zscore]=calRandZscoreForMets(zMetScore,cluster_size,sample)
    %% rpval is pvalues for each enzyme
    % cluster size
    % how many random samples need to be taken
    %%
%     mu=0;
%     sigma=1;
    muk=NaN(numel(cluster_size),1);
    sigmak=NaN(numel(cluster_size),1);
    matrix_zscore=NaN(sample,numel(cluster_size));
    for i=1:numel(cluster_size)
        num=cluster_size(i);
        z_score=zeros(sample,1);
        for k=1:sample
            index=randperm(numel(zMetScore),num);
            met_score=zMetScore(index);
            met_score=met_score(~isnan(met_score));
            
            z_score(k)=sum(met_score)*(1/sqrt(numel(met_score)));
         
        end
        matrix_zscore(:,i)=z_score;   
        muk(i)=mean(z_score);
        sigmak(i)=std(z_score);
    end