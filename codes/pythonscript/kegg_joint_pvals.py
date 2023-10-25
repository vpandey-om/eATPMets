## we want to combine both enrichment analysis
import pandas as pd
import numpy as np
import os
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name
# where the current directory is present.
eatp_metabolism = os.path.dirname(current)
import statsmodels.stats.multitest as mtest ## f
from scipy.stats import combine_pvalues
### combine m9 enrichment based on genes and metabolites
def calcJointPval(metDict,df_gene,df_met,outfile):

    res_df=df_gene.merge(df_met, on='Pathways',suffixes=['_gene','_metabolite'],how='outer')
    res_df['Pvals_gene']=res_df['Pvals_gene'].fillna(1)
    res_df['Pvals_metabolite']=res_df['Pvals_metabolite'].fillna(1)

    # res_df=df_gene.merge(df_met,on='Pathways',suffixes=['_gene','_metabolite'])
    res_df['jointPval']=res_df['Pvals_gene']*res_df['Pvals_metabolite']




    fdr=mtest.multipletests(res_df['jointPval'].values, alpha=0.05, method='fdr_bh')
    res_df['joint_FDR']=fdr[1]
    metlist=res_df['Compounds'].to_list()
    metNamelist=[]

    for item in metlist:
        if type(item) == str:
            tmp=item.split(',')

            aa=[metDict[m.strip()] if m.strip() in metDict.keys() else m.strip() for m in tmp]

            metNamelist.append('|'.join(aa))
        else:
            metNamelist.append('')

    res_df['Compound Names']=metNamelist
    ### filter
    final_df=res_df.sort_values(by=['jointPval'])
    xx=final_df['Number of genes found in pathway']>0
    xx=xx.fillna(0)
    yy=final_df['Number of compounds found in pathway']>0
    yy=yy.fillna(0)
    final_df1=final_df[xx|yy].copy()
    final_df1.to_csv(outfile,sep='\t')


def calcJointPvalGene(df_gene1,df_gene2,outfile):

    res_df=df_gene1.merge(df_gene2, on='Pathways',suffixes=['_1','_2'],how='outer')
    res_df['Pvals_1']=res_df['Pvals_1'].fillna(1)
    res_df['Pvals_2']=res_df['Pvals_2'].fillna(1)
    # res_df=df_gene.merge(df_met,on='Pathways',suffixes=['_gene','_metabolite'])
    res_df['jointPval']=res_df['Pvals_1']*res_df['Pvals_2']
    fdr=mtest.multipletests(res_df['jointPval'].values, alpha=0.05, method='fdr_bh')
    res_df['joint_FDR']=fdr[1]
    final_df=res_df.sort_values(by=['jointPval'])
    ## combine pvals
    jointpval1=[]
    for id in final_df.index:
        tmp=combine_pvalues([final_df.loc[id,'Pvals_1'],final_df.loc[id,'Pvals_2']], method='fisher', weights=None)
        jointpval1.append(tmp[1])
    final_df['jointPval2']=jointpval1
    xx=final_df['Number of genes found in pathway_1']>0
    xx=xx.fillna(False)
    yy=final_df['Number of genes found in pathway_2']>0
    yy=yy.fillna(False)
    final_df1=final_df[xx&yy].copy()
    final_df1.to_csv(outfile,sep='\t')


# in1=os.path.join(eatp_metabolism,'results','M9', 'enrichment_09ST2018_t_180_background_test.txt')
# df_gene1=pd.read_csv(in1,sep='\t')
# in1=os.path.join(eatp_metabolism,'results','M9', 'enrichment_09ST2018_t_180_background_fold_1.5.txt')
# df_gene2=pd.read_csv(in1,sep='\t')
# outfile=os.path.join(eatp_metabolism,'results','M9', 'enrichment_09ST2018_t_180_testJoint.txt')
#
# calcJointPvalGene(df_gene1,df_gene2,outfile)


## first joint pvals M9

df_met1=pd.read_excel(os.path.join(eatp_metabolism,'data', '500uM_ATP_180min_final.xlsx'))

df_met1['KEGG']=df_met1['KEGG'].str.strip()

df_gene=pd.read_csv(os.path.join(eatp_metabolism,'results','M9', 'enrichment_09ST2018_t_180_fold_1.5.txt'),sep='\t')

df_met=pd.read_csv(os.path.join(eatp_metabolism,'results', 'MetaboliteBased', 'enrichment_metabolites_Volcano_M9_1mMATP_pellet.txt'),sep='\t')
metDict=dict(zip(df_met1['KEGG'],df_met1['compound']))


outfile=os.path.join(eatp_metabolism,'results','M9', 'enrichment_combine_09ST2018_t_180_fold_1.5.txt')
calcJointPval(metDict,df_gene,df_met,outfile)
# first joint pvals rich medium
df_met1=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data','Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_1mMATP_180min_vikash.xlsx'))
df_met1['KEGG']=df_met1['KEGG'].str.strip()

df_gene=pd.read_csv(os.path.join(eatp_metabolism,'results','rich', 'enrichment_04ST2018_t_180_fold_1.5.txt'),sep='\t')
df_met=pd.read_csv(os.path.join(eatp_metabolism,'results', 'MetaboliteBased', 'enrichment_metabolites_rich_1mM_ATP_180min.txt'),sep='\t')
metDict=dict(zip(df_met1['KEGG'],df_met1['Column1']))
outfile=os.path.join(eatp_metabolism,'results','M9', 'enrichment_combine_04ST2018_t_180_fold_1.5.txt')
calcJointPval(metDict,df_gene,df_met,outfile)


## joint analysis on rich medium
