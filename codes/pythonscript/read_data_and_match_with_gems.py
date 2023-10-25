## read metabolomics data
import pandas as pd
import cobra
import re

from libsbml import *
reader = SBMLReader()
document = reader.readSBML("/Users/vpandey/projects/gitlabs/eatp_metabolism/models/iML1515.xml")
print(document.getNumErrors())

listOfspecies=document.model.species

mets=[]
metNames=[]
metanetx_ids_all=[]
kegg_ids_all=[]

mets2=[]
metNames2=[]
kegg_ids_all2=[]

for species in listOfspecies:
    if species.annotation_string:
        keggids=re.findall(r'kegg.compound/[CG][0-9]*', species.annotation_string)
        if len(keggids)>0:
            kegg_ids=[ item.replace('kegg.compound/','') for item in keggids ]
            kegg_ids_string='|'.join(kegg_ids)
            for item in kegg_ids:
                kegg_ids_all2.append(item)
                mets2.append(species.id)
                metNames2.append(species.name)

        else:
            kegg_ids_string='NA'
        kegg_ids_all.append(kegg_ids_string)
        metanetids=re.findall(r'metanetx.chemical/MNXM[0-9]*', species.annotation_string)
        if len(metanetids)>0:
            metanet_ids=[ item.replace('metanetx.chemical/','') for item in metanetids ]
            metanet_ids_string='|'.join(metanet_ids)
        else:
            metanet_ids_string='NA'
        metanetx_ids_all.append(metanet_ids_string)
        mets.append(species.id)
        metNames.append(species.name)

df=pd.DataFrame()
df['mets']=mets
df['metNames']=metNames
df['KeggIds']=kegg_ids_all
df['metaNetXIds']=metanetx_ids_all

df2=pd.DataFrame()
df2['mets']=mets2
df2['metNames']=metNames2
df2['KeggIds']=kegg_ids_all2





data100_60=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/100uM_ATP_60min_final.xlsx')
data100_180=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/100uM_ATP_180min_final.xlsx')
data500_60=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/500uM_ATP_60min_final.xlsx')
data500_180=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/500uM_ATP_180min_final.xlsx')
bigg_mets=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/bigg_models_metabolites.xlsx')
df_all=data100_60.rename(columns={"compound": "Compound",'FC':'FC(100uM_60min)','log2(FC)':'log2(FC_100uM_60min)', 'raw.pval':'raw.pval(100uM_60min)','KEGG':'KEGG(100uM_60min)','-log10(p)':'-log10(p_100uM_60min)'})
data500_60=data500_60.rename(columns={"compound": "Compound",'FC':'FC(500uM_60min)','log2(FC)':'log2(FC_500uM_60min)', 'raw.pval':'raw.pval(500uM_60min)','KEGG':'KEGG(500uM_60min)','-log10(p)':'-log10(p_500uM_60min)'})
data100_180=data100_180.rename(columns={"compound": "Compound",'FC':'FC(100uM_180min)','log2(FC)':'log2(FC_100uM_180min)', 'raw.pval':'raw.pval(100uM_180min)','KEGG':'KEGG(100uM_180min)','-log10(p)':'-log10(p_100uM_180min)'})
data500_180=data500_180.rename(columns={"compound": "Compound",'FC':'FC(500uM_180min)','log2(FC)':'log2(FC_500uM_180min)', 'raw.pval':'raw.pval(500uM_180min)','KEGG':'KEGG(500uM_180min)','-log10(p)':'-log10(p_500uM_180min)'})



# df_all=df_all.merge(data100_180, on='compound',suffixes=('(100uM_60min)', '(100uM_180min)'))
# df_all=df_all.merge(data500_60, on='compound',suffixes=('', '(500uM_60min)'))
# df_all=df_all.merge(data500_180, on='compound',suffixes=('', '(500uM_180min)'))


df_all=df_all.merge(data100_180, on='Compound')
df_all=df_all.merge(data500_60, on='Compound')
df_all=df_all.merge(data500_180, on='Compound')

df_all.to_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/all_condition_combined_metabolomics.xlsx')



## get kegg ids form bigg models
bigg_kegg_ids=[]
for link in bigg_mets['database_links'].to_list():
    if type(link) == str:
        keggids=re.findall(r'kegg.compound/[CG][0-9]*', link)

        if len(keggids)>0:
            kegg_ids=[ item.replace('kegg.compound/','') for item in keggids ]

            for item in kegg_ids:
                bigg_kegg_ids.append(item)



##

exp_kegg_ids=[]
names=[]
for i,item, in enumerate(data100_60['KEGG'].to_list()):
    if type(item) == str:
        keggid=re.findall(r'C[0-9]*', item)
        if len(keggid)>0:
            exp_kegg_ids.append(keggid[0])
            names.append(data100_60.loc[i,'Compound'])
df_exp=pd.DataFrame()
df_exp['metNames']=names
df_exp['keggids']=exp_kegg_ids

### found metabolite in the model
found_metabolite=df2[df2['KeggIds'].isin(exp_kegg_ids)]
df_all['KeggIds']=df_all['KEGG(100uM_60min)']
df_all['KeggIds']=[itm[0] if not pd.isnull(itm) else 'NA' for itm in df_all['KeggIds'].str.findall(r'C[0-9]*') ]
###

not_found_in_the_model=set(exp_kegg_ids)-set(kegg_ids_all2)

df_not_found=df_exp[df_exp['keggids'].isin(not_found_in_the_model)]

not_found_all=df_all[df_all['Compound'].isin(df_not_found['metNames'].to_list())]

not_found_all.to_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/not_found_in_model_all.xlsx')


# is there  in bigg
not_found_in_the_model2=not_found_in_the_model-set(bigg_kegg_ids)

df_not_found2=df_exp[df_exp['keggids'].isin(not_found_in_the_model2)]

df_not_found2.to_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/not_found2_in_model.txt',sep='\t')

not_found_all2=df_all[df_all['Compound'].isin(df_not_found2['metNames'].to_list())]

not_found_all2.to_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/not_found_in_model_all2.xlsx')

found_df=pd.merge(found_metabolite,df_all,on='KeggIds',how='left')
found_df.to_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/found_Metabolites_in_model.xlsx')
import pdb; pdb.set_trace()

model=cobra.io.load_matlab_model('/Users/vpandey/projects/gitlabs/eatp_metabolism/models/iML1515.mat')
solution = model.optimize()



## get external ids from the model
