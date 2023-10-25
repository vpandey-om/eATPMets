## we want to combine both enrichment analysis
import pandas as pd
import numpy as np
df1=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/results/metabolite_enrichment.txt',sep='\t')
df2=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/results/rxns_enrichment.txt',sep='\t')
met_df=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/model_mets.xlsx')
metDict=dict(zip(met_df.mets,met_df.metNames))
res_df=df1.merge(df2, on='Pathway',suffixes=['_metabolites','_rxns_genes'])
res_df['jointPval']=res_df[' pvalue_metabolites']*res_df[' pvalue_rxns_genes']
metlist=res_df['DEMetList'].to_list()
metNamelist=[]
for item in metlist:
    if type(item) == str:
        tmp=item.split(',')
        aa=[metDict[m] for m in tmp]
        metNamelist.append('|'.join(aa))
    else:
        metNamelist.append('')

res_df['DE_metabolites']=metNamelist

final_df=res_df.sort_values(by=['jointPval'])

final_df.to_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/results/combined_enrichment_analysis.xlsx',index=None)
