import pandas as pd
import sys,os
from scipy.stats import hypergeom ## for hypergeomtric distributions
import numpy as np
import statsmodels.stats.multitest as mtest ## f

# setting path
# getting the name of the directory
# where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name
# where the current directory is present.
eatp_metabolism = os.path.dirname(current)


class KeggAna:

    def __init__(self):
        self.root_folder=eatp_metabolism

        idToPathDict,pathIdTogeneDict,pathIdToCompound=self.read_eco()
        self.pathwayInfoEco=[idToPathDict,pathIdTogeneDict,pathIdToCompound]
        idToPathDict,pathIdTogeneDict,pathIdToCompound=self.read_ecomodule()
        self.moduleInfoEco=[idToPathDict,pathIdTogeneDict,pathIdToCompound]

    def read_ecomodule(self):
        modules=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'module_name.txt'),sep='\t',header=None)
        geneToModules=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'eco_module.txt'),sep='\t',header=None)

        path_cpd=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'module_compounds.txt'),sep='\t',header=None)

        ##
        idToPathDict=dict(zip(modules[0], modules[1]))
        pathIdTogeneDict={}
        for k,Idx in geneToModules.groupby([1]).indices.items():
            pathname=k.replace('eco_','')
            pathIdTogeneDict[pathname]=[item.replace('eco:','') for item in geneToModules.loc[Idx,0].to_list()]

        ### need to know compunds per paths
        pathIdToCompound={}
        for k,Idx in path_cpd.groupby([0]).indices.items():
            if k in pathIdTogeneDict.keys():
                pathIdToCompound[k]=[item.replace('cpd:','') for item in path_cpd.loc[Idx,1].to_list()]
        # ## get compounds only organism specific
        # # 0) pathway to gene
        # # 1) get gene to ec
        # ## 2) get ec to rxn
        # ## 3) get  rxn to compound
        # for k,Idx in geneToPaths.groupby([1]).indices.items():
        #     geneList=geneToPaths.loc[Idx,0].to_list()
        #     import pdb; pdb.set_trace()

        return idToPathDict,pathIdTogeneDict,pathIdToCompound



    def read_eco(self):
        pathways=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'eco_pathway_list.txt'),sep='\t',header=None)
        geneToPaths=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'pathway_genes.txt'),sep='\t',header=None)
        rxn_enzyme=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'rxn_enzyme.txt'),sep='\t',header=None)
        ec_genes=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'ec_genes.txt'),sep='\t',header=None)
        rxn_cpd=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'reaction_compound.txt'),sep='\t',header=None)
        path_cpd=pd.read_csv(os.path.join(self.root_folder,'Metabolomic_data','kegg', 'pathway_compound.txt'),sep='\t',header=None)

        ##
        idToPathDict=dict(zip(pathways[0], pathways[1]))
        pathIdTogeneDict={}
        for k,Idx in geneToPaths.groupby([1]).indices.items():
            pathIdTogeneDict[k]=[item.replace('eco:','') for item in geneToPaths.loc[Idx,0].to_list()]

        ### need to know compunds per paths
        pathIdToCompound={}
        for k,Idx in path_cpd.groupby([0]).indices.items():
            pathname=k.replace('map','eco')
            if pathname in pathIdTogeneDict.keys():
                pathIdToCompound[pathname]=[item.replace('cpd:','') for item in path_cpd.loc[Idx,1].to_list()]
        # ## get compounds only organism specific
        # # 0) pathway to gene
        # # 1) get gene to ec
        # ## 2) get ec to rxn
        # ## 3) get  rxn to compound
        # for k,Idx in geneToPaths.groupby([1]).indices.items():
        #     geneList=geneToPaths.loc[Idx,0].to_list()
        #     import pdb; pdb.set_trace()

        return idToPathDict,pathIdTogeneDict,pathIdToCompound

    def geneEnrich(self,outfile,targetList,backGroundList,totalGeneNum,method='path'):
        ''' This is the analysis based on genes'''
        if method=='path':
            pathway_to_genes=self.pathwayInfoEco[1]
            pathway_to_Name=self.pathwayInfoEco[0]
        else:
            pathway_to_genes=self.moduleInfoEco[1]
            pathway_to_Name=self.moduleInfoEco[0]




        pvals=[] # list for storing pvalues
        n_list=[] # store number of genes in pathway
        x_list=[] # store number of target genes in pathway
        genes=[] ## we store genes target geens
        pathways=[] ## store pathway
        n_origList=[] ###


        if backGroundList==None:
            M= totalGeneNum
        else:
            M=len(backGroundList)

        ## group clusters
        for k,v in pathway_to_genes.items():
            ## print(k, len(v))


            pathway_genes=set(v) ## genes in the cluster
            if backGroundList==None:
                pathway_background=list(set(v))
            else:
                pathway_background=list(set(backGroundList)& set(v))

            gene_in_pathway= list(set(targetList)& set(v)) ## how many sex-specific genes are found in cluster
            ## how many are in phenoype
            M=M ## total genes
            N=len(targetList) ## sex-specific genes
            n=len(pathway_background) ## total number of cluster genes

            rv = hypergeom(M, n, N) ## hypergeometric  distributions
            x=len(gene_in_pathway) ## number of sex-specific genes in cluster

            ## for p-value we need to get summation of all probalities greater than > x

            pval= 1-rv.cdf(x)

            pvals.append(pval)
            n_list.append(n)
            x_list.append(x)
            genes.append(', '.join(gene_in_pathway))
            n_origList.append(len(pathway_genes))
            pathways.append(pathway_to_Name[k]) ## store cluster number

        fdr=mtest.multipletests(pvals, alpha=0.05, method='fdr_bh')
        enrich={'Pathways':pathways, 'Pvals':pvals,'Pathway size':n_origList,'Pathway size in background':n_list, 'Number of genes found in pathway':x_list,'FDR':fdr[1],'Genes':genes}
        ## create enrichment table
        enrich_df=pd.DataFrame.from_dict(enrich)
        ## filter results
        # enrich_df=enrich_df[enrich_df['Number of genes found in pathway']>=1]
        enrich_soreted=enrich_df.sort_values(by=['Pvals'],ascending=True)
        enrich_soreted.to_csv(outfile,sep='\t',index=None)


    def cpdEnrich(self,outfile,targetList,backGroundList,totalMetNum):
        ''' This is the analysis based on genes'''

        pathway_to_cpds=self.pathwayInfoEco[2]
        pathway_to_Name=self.pathwayInfoEco[0]

        pvals=[] # list for storing pvalues
        n_list=[] # store number of genes in pathway
        x_list=[] # store number of target genes in pathway
        cpds=[] ## we store genes target geens
        pathways=[] ## store pathway
        n_origList=[] ###


        if backGroundList==None:
            M= totalMetNum
        else:
            M=len(backGroundList)

        ## group clusters
        for k,v in pathway_to_cpds.items():
            ## print(k, len(v))
            # pathways.append(pathway_to_Name[k]) ## store cluster number

            pathway_cpds=set(v) ## genes in the cluster
            if backGroundList==None:
                pathway_background=list(set(v))
            else:
                pathway_background=list(set(backGroundList)& set(v))

            cpd_in_pathway= list(set(targetList)& set(v)) ## how many sex-specific genes are found in cluster
            ## how many are in phenoype
            M=M ## total genes
            N=len(targetList) ## sex-specific genes
            n=len(pathway_background) ## total number of cluster genes

            rv = hypergeom(M, n, N) ## hypergeometric  distributions
            x=len(cpd_in_pathway) ## number of sex-specific genes in cluster

            ## for p-value we need to get summation of all probalities greater than > x

            pval= 1-rv.cdf(x)
            if np.isnan(pval):
                import pdb; pdb.set_trace()
                print(pval)
            else:
                pvals.append(pval)
                n_list.append(n)
                x_list.append(x)
                cpds.append(', '.join(cpd_in_pathway))
                n_origList.append(len(pathway_cpds))
                pathways.append(pathway_to_Name[k]) ## store cluster number
        fdr=mtest.multipletests(pvals, alpha=0.05, method='fdr_bh')
        enrich={'Pathways':pathways, 'Pvals':pvals,'Pathway size':n_origList,'Pathway size in background':n_list, 'Number of compounds found in pathway':x_list,'FDR':fdr[1],'Compounds':cpds}
        ## create enrichment table
        enrich_df=pd.DataFrame.from_dict(enrich)
        ## filter results
        # enrich_df=enrich_df[enrich_df['Number of compounds found in pathway']>=1]
        enrich_soreted=enrich_df.sort_values(by=['Pvals'],ascending=True)
        enrich_soreted.to_csv(outfile,sep='\t',index=None)



## read M9 genetic data

kd=KeggAna()
# #
df=pd.read_excel(os.path.join(kd.root_folder,'data', '09ST2018_Summary.xlsx'),sheet_name='Vikash_logFC')
# df=pd.read_csv(os.path.join(kd.root_folder,'data', '09ST2018_Summary.txt'),sep='\t')
df.to_excel(os.path.join(kd.root_folder,'data', '09ST2018_Summary_reduced.xlsx'))
logfold=np.log2(3/2)
filter_df=df[abs(df['log2_t180'])>=logfold]

filter_df.to_csv(os.path.join(kd.root_folder,'data', '09ST2018_filter.txt'),sep='\t')
deregulated_genes=filter_df['Uniprot'].dropna().to_list()
background_genes=df['Uniprot'].dropna().to_list()
print(len(deregulated_genes),len(background_genes))
outfile=os.path.join(kd.root_folder,'results', 'M9','enrichment_09ST2018_t_180_background_fold_1.5.txt')

kd.geneEnrich(outfile,deregulated_genes,background_genes,len(background_genes))

outfile=os.path.join(kd.root_folder,'results', 'M9','enrichment_09ST2018_t_180_background_test.txt')

kd.geneEnrich(outfile,background_genes,None,4736)

outfile=os.path.join(kd.root_folder,'results', 'M9', 'enrichment_09ST2018_t_180_fold_1.5.txt')
kd.geneEnrich(outfile,deregulated_genes,None,4736)
# outfile=os.path.join(kd.root_folder,'results', 'enrichment_09ST2018_t_180_module.txt')
# kd.geneEnrich(outfile,deregulated_genes,None,4736,'module')




df=pd.read_excel(os.path.join(kd.root_folder,'data', 'All_gfp_vfinal_04ST2018.xlsx'),sheet_name='Vikash_logFC')
# df=pd.read_csv(os.path.join(kd.root_folder,'data', '09ST2018_Summary.txt'),sep='\t')
df.to_excel(os.path.join(kd.root_folder,'data', 'All_gfp_vfinal_04ST2018_reduced.xlsx'))


logfold=np.log2(3/2)
filter_df=df[abs(df['log2_t180'])>=logfold]
filter_df.to_csv(os.path.join(kd.root_folder,'data', '04ST2018_filter.txt'),sep='\t')
deregulated_genes=filter_df['Gene_number'].dropna().to_list()
background_genes=df['Gene_number'].dropna().to_list()

print(len(deregulated_genes),len(background_genes))
outfile=os.path.join(kd.root_folder,'results', 'rich','enrichment_04ST2018_t_180_background_fold_1.5.txt')


kd.geneEnrich(outfile,deregulated_genes,background_genes,len(background_genes))
outfile=os.path.join(kd.root_folder,'results', 'rich', 'enrichment_04ST2018_t_180_fold_1.5.txt')
kd.geneEnrich(outfile,deregulated_genes,None,4736)
# outfile=os.path.join(kd.root_folder,'results', 'enrichment_09ST2018_t_180_module.txt')
# kd.geneEnrich(outfile,deregulated_genes,None,4736,'module')




### based on metabolites

## read metabolite data
logfold=np.log2(3/2)
df_met=pd.read_excel(os.path.join(kd.root_folder,'data', '100uM_ATP_180min_final.xlsx'))
df_met['KEGG']=df_met['KEGG'].str.strip()

filter_df=df_met[abs(df_met['log2(FC)'])>logfold]
deregulated_mets=filter_df['KEGG'].dropna().to_list()
background_mets=df_met['KEGG'].dropna().to_list()

print(len(deregulated_mets),len(background_mets))
outfile=os.path.join(kd.root_folder,'results', 'MetaboliteBased', 'enrichment_metabolites_M9_100uM_ATP_180min.txt')
# kd.cpdEnrich(outfile,deregulated_mets,background_mets,len(background_mets))
kd.cpdEnrich(outfile,deregulated_mets,None,3000)


logfold=np.log2(3/2)
df_met1=pd.read_excel(os.path.join(kd.root_folder,'data', '500uM_ATP_180min_final.xlsx'))

df_met1['KEGG']=df_met1['KEGG'].str.strip()

filter_df=df_met1[abs(df_met1['log2(FC)'])>logfold]
deregulated_mets=filter_df['KEGG'].dropna().to_list()
background_mets=df_met1['KEGG'].dropna().to_list()

print(len(deregulated_mets),len(background_mets))
outfile=os.path.join(kd.root_folder,'results', 'MetaboliteBased', 'enrichment_metabolites_M9_500uM_ATP_180min.txt')
# kd.cpdEnrich(outfile,deregulated_mets,background_mets,len(background_mets))
kd.cpdEnrich(outfile,deregulated_mets,None,3000)

#
# logfold=np.log2(3/2)
#
#
# df_met1=pd.read_excel(os.path.join(kd.root_folder,'Metabolomic_data', 'M9_pellet_1mMATP_60min_180min_vikash.xlsx'))
#
# df_met1['KEGG']=df_met1['KEGG'].str.strip()
#
# filter_df=df_met1[abs(df_met1['log2(FC)'])>logfold]
# filter_df.to_csv(os.path.join(kd.root_folder,'data', 'M9_1mMATP_60min_180min_filter.txt'),sep='\t')
# deregulated_mets=filter_df['KEGG'].dropna().to_list()
# background_mets=df_met1['KEGG'].dropna().to_list()
#
# print(len(deregulated_mets),len(background_mets))
# outfile=os.path.join(kd.root_folder,'results', 'MetaboliteBased', 'enrichment_metabolites_1mMATP_60min_180min.txt')
# # kd.cpdEnrich(outfile,deregulated_mets,background_mets,len(background_mets))
# kd.cpdEnrich(outfile,deregulated_mets,None,3000)



logfold=np.log2(3/2)


df_met1=pd.read_excel(os.path.join(kd.root_folder,'Metabolomic_data', 'Volcano_M9_1mMATP_pellet_1.xlsx'))
df_met1['KEGG']=df_met1['KEGG'].str.strip()

filter_df=df_met1[abs(df_met1['log2(FC)'])>logfold]



filter_df.to_csv(os.path.join(kd.root_folder,'data', 'Volcano_M9_1mMATP_pellet.txt'),sep='\t')
deregulated_mets=filter_df['KEGG'].dropna().to_list()
background_mets=df_met1['KEGG'].dropna().to_list()
print(len(deregulated_mets),len(background_mets))
outfile=os.path.join(kd.root_folder,'results', 'MetaboliteBased', 'enrichment_metabolites_Volcano_M9_1mMATP_pellet.txt')
# kd.cpdEnrich(outfile,deregulated_mets,background_mets,len(background_mets))
kd.cpdEnrich(outfile,deregulated_mets,None,3000)




### Rich medium metabolites

logfold=np.log2(3/2)
df_met1=pd.read_excel(os.path.join(kd.root_folder,'Metabolomic_data', 'Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_1mMATP_180min_vikash.xlsx'))

df_met1['KEGG']=df_met1['KEGG'].str.strip()

filter_df=df_met1[abs(df_met1['log2(FC)'])>logfold]
filter_df.to_csv(os.path.join(kd.root_folder,'data', 'rich_1mM_180min_filter.txt'),sep ='\t')
deregulated_mets=filter_df['KEGG'].dropna().to_list()
background_mets=df_met1['KEGG'].dropna().to_list()

print(len(deregulated_mets),len(background_mets))
outfile=os.path.join(kd.root_folder,'results', 'MetaboliteBased', 'enrichment_metabolites_rich_1mM_ATP_180min.txt')
# kd.cpdEnrich(outfile,deregulated_mets,background_mets,len(background_mets))
kd.cpdEnrich(outfile,deregulated_mets,None,3000)
print('rich')
print(deregulated_mets)
## Rich medium


logfold=np.log2(3/2)
df_met1=pd.read_excel(os.path.join(kd.root_folder,'Metabolomic_data', 'Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_500uM_ATP_180min_vikash.xlsx'))

df_met1['KEGG']=df_met1['KEGG'].str.strip()

filter_df=df_met1[abs(df_met1['log2(FC)'])>logfold]
deregulated_mets=filter_df['KEGG'].dropna().to_list()
background_mets=df_met1['KEGG'].dropna().to_list()

print(len(deregulated_mets),len(background_mets))
outfile=os.path.join(kd.root_folder,'results', 'MetaboliteBased', 'enrichment_metabolites_rich_500uM_ATP_180min.txt')
# kd.cpdEnrich(outfile,deregulated_mets,background_mets,len(background_mets))
kd.cpdEnrich(outfile,deregulated_mets,None,3000)
