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


### pathway enrichment based on Flux

infile=os.path.join(eatp_metabolism,'results','modelSim', 'M9_rxns_enrichment.txt')
deregCol='DERxns'
pvalCol=' pvalue'
pathwayCol="Subsystems"
outpdf=os.path.join(eatp_metabolism,'results','modelSim', 'M9_rxns_enrichment.pdf')
title='M9 rxn flux based enrichment'
pvalcut=0.05
pathway_size='Reactions'

plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)




infile=os.path.join(eatp_metabolism,'results','modelSim', 'LB_rxns_enrichment.txt')
deregCol='DERxns'
pvalCol=' pvalue'
pathwayCol="Subsystems"
outpdf=os.path.join(eatp_metabolism,'results','modelSim', 'LB_rxns_enrichment.pdf')
title='LB rxn flux based enrichment'
pvalcut=0.05
pathway_size='Reactions'

plot_enrichmentFile(infile, deregCol,pvalCol, pathwayCol,outpdf,pvalcut,title,pathway_size)
