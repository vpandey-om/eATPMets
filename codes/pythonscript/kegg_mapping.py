import pandas as pd
met_df=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/100uM_ATP_180min_final.xlsx')
gene_df=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/M9_zazlaver_logFC.xlsx',sheet_name='logFC_expressed')
cut=1.5
dereg_met=met_df[(met_df['FC']>cut) | (met_df['FC']<(1/cut))]
keggids=dereg_met['KEGG'].str.strip().to_list()
dereg_genes_df=gene_df[(gene_df['log2_t180']>cut) | (gene_df['log2_t180']>(1/cut))]
genes=dereg_genes_df['Gene_name3'].to_list()
all_genes=[]
for g in genes:
    if type(g) == str:
        aa=g.split(',')
        for a in aa:
            all_genes.append('eco:'+a.strip())
## get gene ids eco:b0114
# make for kegg search
all_entries=all_genes+keggids
df=pd.DataFrame()
df['allIds']=all_entries
df.to_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/data/deregulated_ids.txt',index=None,header=None)
