import pandas as pd ### based on kegg compounds
import numpy as np
import os
import pickle

logfold=np.log2(3/2)

# setting path
# getting the name of the directory
# where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name
# where the current directory is present.
eatp_metabolism = os.path.dirname(current)

listofdata=[]
listofexp=[]
cut=1.5
cut2=0.06
logfold=np.log2(cut)
df_met1=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data', 'Volcano_M9_1mMATP_pellet_1.xlsx'))
df_met1['KEGG']=df_met1['KEGG'].str.strip()

listofdata.append(df_met1)
listofexp.append('M9_metabolite')


# m9_met_df=df_met1[abs(df_met1['log2(FC)'])>logfold]
# m9_met_df_p=df_met1[df_met1['raw,pval']<0.06]

# m9_met_df=pd.concat([m9_met_df,m9_met_df_p],axis=0).drop_duplicates()
m9_met_df=df_met1[(abs(df_met1['log2(FC)'])>logfold) | (df_met1['raw,pval']<cut2)]
m9_met_df.to_excel(os.path.join(eatp_metabolism,'Metabolomic_data', 'M9_1mMATP_deregulated.xlsx'))

logfold=np.log2(cut)
df_met1=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data', 'Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_1mMATP_180min_vikash.xlsx'))
df_met1['KEGG']=df_met1['KEGG'].str.strip()

listofdata.append(df_met1)
listofexp.append('rich_metabolite')


# rich_met_df=df_met1[abs(df_met1['log2(FC)'])>logfold]
# rich_met_df_p=df_met1[df_met1['raw.pval']<0.06]
# rich_met_df=pd.concat([rich_met_df,rich_met_df_p],axis=0).drop_duplicates()
rich_met_df=df_met1[(abs(df_met1['log2(FC)'])>logfold) | (df_met1['raw.pval']<cut2)]
rich_met_df.to_excel(os.path.join(eatp_metabolism,'Metabolomic_data', 'rich_1mMATP_deregulated.xlsx'))



# rich_met_df=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/rich_1mM_180min_filter.txt',sep='\t')
rich_gene_df=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/04ST2018_filter.txt',sep='\t')
listofdata.append(rich_gene_df)
listofexp.append('rich_gene')

# m9_met_df=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/Volcano_M9_1mMATP_pellet.txt',sep='\t')
m9_gene_df=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/09ST2018_filter.txt',sep='\t')

listofdata.append(m9_gene_df)
listofexp.append('m9_genes')

# pickle.dump([listofdata,listofexp],open('/Users/vikash/gitlab/pathcompnet/data/sophie_data.pickle','wb'))


M9mets=m9_met_df['KEGG'].dropna().to_list()
richmets=rich_met_df['KEGG'].dropna().to_list()
M9=['C00575', 'C15517', 'C00147', 'C00739', 'C08277', 'C00212',
'C00127', 'C05512', 'C00117', 'C21016', 'C00624', 'C01042', 'C00105', 'C00294']
rich=['C15517', 'C00294', 'C00328', 'C00078', 'C21016', 'C00029']

## all metabolites

# for item in M9:
#     print(item+'\tred')
rich_color_df=pd.DataFrame()

rich_color_df['items']=rich_gene_df['Gene_number'].to_list()+richmets
rich_color_df['color']='red'
rich_color_df.to_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/colorKEGG/rich_kegg_color.txt',sep='\t',index=None,header=None)
m9_color_df=pd.DataFrame()
m9_color_df['items']=m9_gene_df['Uniprot'].to_list()+M9mets
m9_color_df['color']='red'
m9_color_df.to_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/colorKEGG/m9_kegg_color.txt',sep='\t',index=None,header=None)

###
models=pd.read_excel('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/found_Metabolites_in_model.xlsx')

#### M9 model whole data



####

keggids=set(models['KeggIds'].to_list())

commonM9=set(M9mets) & keggids
commonrich=set(richmets) & keggids
df=models[models['KeggIds'].isin(commonM9)]

df1=df[['mets','metNames','KeggIds']]
df1=df1.rename(columns={'KeggIds':'KEGG'})
final_m9_reg_df=df1.merge(m9_met_df,on='KEGG')
final_m9_reg_df.to_excel('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/M9_dereg_model_mets.xlsx')
df=models[models['KeggIds'].isin(commonrich)]

df1=df[['mets','metNames','KeggIds']]
df1=df1.rename(columns={'KeggIds':'KEGG'})
final_rich_reg_df=df1.merge(rich_met_df,on='KEGG')
final_rich_reg_df.to_excel('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/rich_dereg_model_mets.xlsx')
# rich_gene_df.to_excel('/Users/vikash/Documents/MATLAB/eatp_metabolism/data/rich_dereg_model_genes.xlsx')


import pdb; pdb.set_trace()
print(df['metNames'].unique())


print(df['metNames'].unique())

import pdb; pdb.set_trace()
