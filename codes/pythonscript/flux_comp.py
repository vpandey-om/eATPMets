import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# flux_df=pd.read_excel('/Users/vpandey/projects/gitlabs/eatp_metabolism/results/flux_result.xlsx',header=None)
flux=[]
type=[]
rxns=[]

df=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/results/LB_REMI.txt',header=None,sep='\t')
rxn_pathway=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/matfiles/rxns_flux_enrichment.txt',sep='\t')
rxns=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/matfiles/lb_dere_rxns.txt',sep='\t',header=None)
subys_rxn=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/matfiles/subsysm_rxns.txt',sep='\t',header=None)
import pdb; pdb.set_trace()
flux_data=df.values

xx=['PPM2','IMPD','DRPA','GLCt2pp','SULabcpp']
## filter flux df
data=flux_df.iloc[:,1:61].copy()


order_rxns=xx

for ind in flux_df.index:

    if flux_df.iloc[ind,0] in xx:

        for item in flux_df.iloc[ind,1:31]:
            flux.append(item)
            rxns.append(flux_df.iloc[ind,0])
            type.append('WT')
        for item in flux_df.iloc[ind,31:61]:
            flux.append(item)
            rxns.append(flux_df.iloc[ind,0])
            type.append('eATP')

resdf=pd.DataFrame()
resdf['flux']=flux
resdf['type']=type
resdf['Reactions']=rxns


plt.figure()

sns.violinplot(x="Reactions", y="flux", data=resdf,inner='box',
                     hue='type', order=order_rxns) ## 'quartile' 'box'
# sns.swarmplot(y="Relative oocyst growth",
#                 x="oo",data=oocys_df,edgecolor="gray",size=4,color='#FC8D62',order=['Female','Male','Both sexes'])
# sns.stripplot(y="flux",
#                 x="Reactions",hue='type',data=resdf,edgecolor="gray",size=8,order=order_rxns)

#sns.stripplot(x="Gametocyte", y="Relative female fertility", data=female_tmp,palette=color_gam_dict, hue='Phenotypes',edgecolor="gray",size=8)


# figure = sns_plot.get_figure()
plt.tick_params(axis='x', rotation=20)
# plt.xlabel('Gametocyte', fontsize=20,labelpad=10)

plt.savefig('/Users/vpandey/projects/gitlabs/eatp_metabolism/results/flux_rxns.pdf',
            format='pdf',dpi=300)

import pdb; pdb.set_trace()
plt.close()
