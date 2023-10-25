import pandas as pd
import sys,os
import numpy as np
import cobra
import re

# setting path
# getting the name of the directory
# where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name
# where the current directory is present.
eatp_metabolism = os.path.dirname(current)



from libsbml import *
reader = SBMLReader()
document = reader.readSBML(os.path.join(eatp_metabolism,'models', 'iML1515.xml'))
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

bigg_mets=pd.read_excel(os.path.join(eatp_metabolism,'data', 'bigg_models_metabolites.xlsx'))
M9_mets=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data', 'M9_pellet_1mMATP_60min_180min_vikash.xlsx'))
LB_mets=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data', 'Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_1mMATP_180min_vikash.xlsx'))


## get kegg ids form bigg models
bigg_kegg_ids=[]
for link in bigg_mets['database_links'].to_list():
    if type(link) == str:
        keggids=re.findall(r'kegg.compound/[CG][0-9]*', link)

        if len(keggids)>0:
            kegg_ids=[ item.replace('kegg.compound/','') for item in keggids ]

            for item in kegg_ids:
                bigg_kegg_ids.append(item)

###### For M9 matching
exp_kegg_ids=[]
names=[]
for i,item, in enumerate(M9_mets['KEGG'].to_list()):
    if type(item) == str:
        keggid=re.findall(r'C[0-9]*', item)
        if len(keggid)>0:
            exp_kegg_ids.append(keggid[0])
            names.append(M9_mets.loc[i,'Compound'])
df_exp=pd.DataFrame()
df_exp['metNames']=names
df_exp['keggids']=exp_kegg_ids

### found metabolite in the model
found_metabolite=df2[df2['KeggIds'].isin(exp_kegg_ids)]
M9_mets['KeggIds']=M9_mets['KEGG']
M9_mets['KeggIds']=[itm[0] if not pd.isnull(itm) else 'NA' for itm in M9_mets['KeggIds'].str.findall(r'C[0-9]*') ]
found_df=pd.merge(found_metabolite,M9_mets,on='KeggIds',how='left')
found_df.to_excel(os.path.join(eatp_metabolism,'data','M9_match_model_metabolite_data.xlsx'),index=False)
print( len(found_df['KeggIds'].unique()))
### for LB data maching

exp_kegg_ids=[]
names=[]
for i,item, in enumerate(LB_mets['KEGG'].to_list()):
    if type(item) == str:
        keggid=re.findall(r'C[0-9]*', item)
        if len(keggid)>0:
            exp_kegg_ids.append(keggid[0])
            names.append(M9_mets.loc[i,'Compound'])
df_exp=pd.DataFrame()
df_exp['metNames']=names
df_exp['keggids']=exp_kegg_ids

### found metabolite in the model
found_metabolite=df2[df2['KeggIds'].isin(exp_kegg_ids)]
LB_mets['KeggIds']=LB_mets['KEGG']
LB_mets['KeggIds']=[itm[0] if not pd.isnull(itm) else 'NA' for itm in LB_mets['KeggIds'].str.findall(r'C[0-9]*') ]
found_df=pd.merge(found_metabolite,LB_mets,on='KeggIds',how='left')
found_df.to_excel(os.path.join(eatp_metabolism,'data','LB_match_model_metabolite_data.xlsx'),index=False)
print( len(found_df['KeggIds'].unique()))
