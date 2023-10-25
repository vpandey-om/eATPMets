import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import *
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name
# where the current directory is present.
eatp_metabolism = os.path.dirname(current)



def plot_enrichmentFile(infile, deregCol,pvalCol,pathwayCol,outpdf='output.pdf',pvalcut=0.05,title='title',pathway_size='Pathway size'):
    ''' here we are plotting enrichment'''
    df_gene=pd.read_csv(infile,sep='\t')
    filter_df=df_gene[df_gene[pathway_size]>1]
    filter_df=filter_df[filter_df[deregCol]>=1]
    filter_df=filter_df[filter_df[pvalCol]<=pvalcut]
    filter_df['richFactor']=filter_df[deregCol]/filter_df[pathway_size]
    # scaleFactor = 1500/filter_df['Pathway size'].max()
    # sizes=(filter_df['Pathway size'].min()*scaleFactor,filter_df['Pathway size'].max()*scaleFactor)
    filter_df[pathwayCol]=filter_df[pathwayCol].str.replace('- Escherichia coli K-12 MG1655','')
    log10pval='-log10'+pvalCol
    filter_df[log10pval]=-np.log10(filter_df[pvalCol])
    # fig = plt.figure(figsize=(30, 8))
    ##### see with plotnine
    plots =(ggplot(aes(x='richFactor', y=pathwayCol), filter_df)+geom_point(aes(size=pathway_size,color=log10pval))
    + scale_color_gradient(low='blue', high='red')+ ggtitle(title))
    plots.save(outpdf)



def plot_enrichmentFileJoint(infile, deregCol1,deregCol2,pvalCol,pathwayCol,pathway_size1,pathway_size2,outpdf='output.pdf',pvalcut=0.05,title='title'):
    ''' here we are plotting enrichment'''
    df_gene=pd.read_csv(infile,sep='\t')
    filter_df=df_gene[df_gene[pathway_size1]>1]
    filter_df=df_gene[df_gene[pathway_size2]>1]
    filter_df=filter_df[filter_df[deregCol1]>=1]
    filter_df=filter_df[filter_df[deregCol2]>=1]
    filter_df=filter_df[filter_df[pvalCol]<=pvalcut]
    filter_df['richFactor1']=filter_df[deregCol1]/filter_df[pathway_size1]
    filter_df['richFactor2']=filter_df[deregCol2]/filter_df[pathway_size2]
    filter_df['richFactor']=filter_df['richFactor1']+filter_df['richFactor2']
    # scaleFactor = 1500/filter_df['Pathway size'].max()
    # sizes=(filter_df['Pathway size'].min()*scaleFactor,filter_df['Pathway size'].max()*scaleFactor)
    filter_df['Pathway size']=filter_df[pathway_size1]
    filter_df[pathwayCol]=filter_df[pathwayCol].str.replace('- Escherichia coli K-12 MG1655','')
    log10pval='-log10'+pvalCol
    filter_df[log10pval]=-np.log10(filter_df[pvalCol])
    # fig = plt.figure(figsize=(30, 8))
    ##### see with plotnine
    plots =(ggplot(aes(x='richFactor', y='Pathways'), filter_df)+geom_point(aes(size='Pathway size',color=log10pval))
    + scale_color_gradient(low='blue', high='red')+ ggtitle(title))
    plots.save(outpdf)


### pathway enrichment based on Flux
infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'rxns_flux_enrichment.txt')
deregCol='DERxns'
pvalCol=' pvalue'
pathwayCol="Pathway"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'enrichment_rxns_flux_enrichment.pdf')
title='Lb rxn flux based enrichment'
pvalcut=0.05
pathway_size='TRxns'

plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)


#### Figure For M9 genetic screen enrichment

infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_09ST2018_t_180_fold_1.5.txt')
deregCol='Number of genes found in pathway'
pvalCol='Pvals'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_09ST2018_t_180_fold_1.5.pdf')
title='M9 gene based enrichment'
pvalcut=0.05
pathway_size='Pathway size'

plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)

#### Figure For rich medium(LB) genetic screen enrichment

infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_04ST2018_t_180_fold_1.5.txt')
deregCol='Number of genes found in pathway'
pvalCol='Pvals'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_04ST2018_t_180_fold_1.5.pdf')
title='rich medium(LB) gene based enrichment'
pvalcut=0.05
pathway_size='Pathway size'
plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)

#### Metabolite based enrichments ######
## for M9 100 and 500 at 180 minutes
# for M9 100uM at 180 minutes


infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_metabolites_M9_100uM_ATP_180min.txt')
deregCol='Number of compounds found in pathway'
pvalCol='FDR'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_metabolites_M9_100uM_ATP_180min.pdf')
title='M9 100uM at 180 minutes metabolite enrichment'
pvalcut=0.05
pathway_size='Pathway size'
plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)


infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_metabolites_M9_500uM_ATP_180min.txt')
deregCol='Number of compounds found in pathway'
pvalCol='FDR'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_metabolites_M9_500uM_ATP_180min.pdf')
title='M9 500uM at 180 minutes metabolite enrichment'
pvalcut=0.05
pathway_size='Pathway size'
plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)


## for LB 1000 and 500 at 180 minutes
# for LB 500uM at 180 minutes

infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_metabolites_rich_500uM_ATP_180min.txt')
deregCol='Number of compounds found in pathway'
pvalCol='FDR'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_metabolites_rich_500uM_ATP_180min.pdf')
title='rich 500uM at 180 minutes metabolite enrichment'
pvalcut=0.05
pathway_size='Pathway size'
plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)


infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_metabolites_rich_1mM_ATP_180min.txt')
deregCol='Number of compounds found in pathway'
pvalCol='FDR'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_metabolites_rich_1mM_ATP_180min.pdf')
title='rich 1mM at 180 minutes metabolite enrichment'
pvalcut=0.05
pathway_size='Pathway size'
plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)

# crazy compare for LB 1mM 60 and 180 minutes

infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_metabolites_Volcano_M9_1mMATP_pellet.txt')
deregCol='Number of compounds found in pathway'
pvalCol='FDR'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_metabolites_Volcano_M9_1mMATP_pellet.pdf')
title='enrichment_metabolites_M9_1mMATP'
pvalcut=0.05
pathway_size='Pathway size'
plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)

### combined analysis for M9 : based on genes and metabolites
infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_combine_09ST2018_t_180_fold_1.5.txt')
deregCol1='Number of compounds found in pathway'
deregCol2='Number of genes found in pathway'
pvalCol='joint_FDR'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_combine_09ST2018_t_180_fold_1.5.pdf')
title='M9 combine gene and metabolite enrichment'
pvalcut=0.05
pathway_size1='Pathway size_metabolite'
pathway_size2='Pathway size_gene'
plot_enrichmentFileJoint(infile, deregCol1,deregCol2,pvalCol,pathwayCol,pathway_size1,pathway_size2,outpdf,pvalcut,title)


infile=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Files', 'enrichment_combine_04ST2018_t_180_fold_1.5.txt')
deregCol1='Number of compounds found in pathway'
deregCol2='Number of genes found in pathway'
pvalCol='joint_FDR'
pathwayCol="Pathways"
outpdf=os.path.join(eatp_metabolism,'results','EnrichmentResult', 'Figures', 'enrichment_combine_04ST2018_t_180_fold_1.5.pdf')
title='rich combine gene and metabolite enrichment'
pvalcut=0.05
pathway_size1='Pathway size_metabolite'
pathway_size2='Pathway size_gene'
plot_enrichmentFileJoint(infile, deregCol1,deregCol2,pvalCol,pathwayCol,pathway_size1,pathway_size2,outpdf,pvalcut,title)





# ax=sns.scatterplot(data=filter_df, x="rich_factor", y="Pathways", size='Pathway size', legend=False, sizes=sizes,hue='log10pval',palette="Reds")
#
# # # show the graph
# # norm = plt.Normalize(filter_df["log10pval"].min(), filter_df["log10pval"].max())
# # sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
# # sm.set_array([])
# #
# # # Remove the legend and add a colorbar
# # ax.get_legend().remove()
# # ax.figure.colorbar(sm)
# ax.set_yticklabels([])
# outpdf=os.path.join(eatp_metabolism,'results', 'M9_genetic_enrichemnt_points.pdf')
# # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
# plt.savefig(outpdf)
