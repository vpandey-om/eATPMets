import pandas as pd
# import umap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans,AffinityPropagation
import numpy as np
from bioinfokit.visuz import cluster
import pandas.api.types as pdtypes
import plotnine as pn

def draw_pca(u,n_components=3):
    fig = plt.figure()
    if n_components == 2:
        ax = fig.add_subplot(111)
        ax.scatter(u[:,0], u[:,1], c=colors,s=3)
        for i in range(u.shape[0]):
            ax.annotate(str(i), (u[i,0], u[i,1]), fontsize=2)
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(u[:,0], u[:,1], u[:,2], c=colors, s=10)
    # plt.title(title, fontsize=18)
    plt.savefig('testing.pdf')






#df=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/results/LB_REMI.txt',header=None,sep='\t')
df=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/LB_REMI.txt',header=None,sep='\t')
data=df.values
#subys_rxn=pd.read_csv('/Users/vpandey/projects/gitlabs/eatp_metabolism/matfiles/subsysm_rxns.txt',sep='\t',header=None)
subys_rxn=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/subsysm_rxns.txt',sep='\t',header=None)
altsol=500


############## plot violin for which fluxes are differed sigficnatly
df_path=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/matfiles/rxns_flux_enrichment.txt',sep='\t')
derxn_subsystem=df_path[df_path['DERxns']>0]

for index in derxn_subsystem.index:
    subsystem=derxn_subsystem.loc[index,'Pathway']
    derxns=derxn_subsystem.loc[index,'DErxnsList'].split(',')
    tmp=subys_rxn[subys_rxn[0].isin(derxns)]
    # get control flux
    # tmp= tmp[tmp[0].str.contains("O2tpp") == False]

    if not tmp.empty:
        ### remove 02tpp transport
        ctlArr=data[tmp.index,0:altsol]
        eatpArr=data[tmp.index,(altsol+1):(2*altsol+1)]
        diffFlux=np.mean(abs(eatpArr-ctlArr),axis=1)
        rxns=tmp.loc[tmp.index,0].to_list()
        ## make dataframe for plots
        plotdf=pd.DataFrame()
        ctlFluxAllSol=ctlArr.flatten()
        eatpFluxAllSol=eatpArr.flatten()
        allRxns=[[item]*altsol for item in rxns]
        flat_allRxns = [item for sublist in allRxns for item in sublist]
        conds=['Control','eATP']
        allconds=[[item]*len(flat_allRxns) for item in conds]
        flat_allconds = [item for sublist in allconds for item in sublist]
        plotdf['Reactions']=flat_allRxns+flat_allRxns
        plotdf['Condition']=flat_allconds
        plotdf['Flux']=np.concatenate((ctlFluxAllSol,eatpFluxAllSol), axis=None)
        plotdf['Condition'] = plotdf['Condition'].astype(pdtypes.CategoricalDtype(categories=conds))
        plotdf['Reactions'] = plotdf['Reactions'].astype(pdtypes.CategoricalDtype(categories=rxns))
        ## plot_for_each_reactions

        for rxn in rxns:
            tmpdf=plotdf[plotdf['Reactions'].isin([rxn])]
            plots=(pn.ggplot(tmpdf, pn.aes(x='Reactions', y='Flux', fill='Condition'))
                # + pn.coord_flip()
            + pn.geom_violin(tmpdf, style='full')
            # + pn.geom_point(color='none', size=2, show_legend=False)

            + pn.scale_fill_manual(values=['dodgerblue', 'darkorange'])
            + pn.theme_classic()
            + pn.labs(title=subsystem)
            )
            subsystem=subsystem.replace('/','_')
            plots.save('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/pca/flux_plot/'+subsystem+'_'+rxn+'.pdf')
import pdb; pdb.set_trace()







diff_flux=data[:,(altsol+1):(2*altsol+1)]-data[:,0:altsol]
## groupby  subsytems
abs_diff_flux=abs(diff_flux)
subsystems=[]
pathway_flux=np.zeros((len(subys_rxn[1].unique()),altsol))

counter=0
for subsys,idx in subys_rxn.groupby([1]).indices.items():
     #rxn and order are same
     pathway_flux[counter,:]=diff_flux[idx,:].mean(axis=0)
     counter=counter+1
     subsystems.append(subsys)





pathway_fluxT=pathway_flux.T
pca = PCA(n_components=30, svd_solver='full')
pca.fit(pathway_fluxT)
pca_score = pca.transform(pathway_fluxT)
print(pca.explained_variance_ratio_)
print(np.cumsum(pca.explained_variance_ratio_))
loadings = pca.components_
num_pc = pca.n_features_
pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
# loadings_df['variable'] = ['alt'+str(i) for i in range(altsol)]
loadings_df['variable'] = subsystems
loadings_df = loadings_df.set_index('variable')

## scree plot

cluster.screeplot(obj=[pc_list[:30], pca.explained_variance_ratio_],figtype='pdf')

## pca plot


cluster.pcaplot(x=loadings[0], y=loadings[1], labels=loadings_df.index,
    var1=round(pca.explained_variance_ratio_[0]*100, 2),
    var2=round(pca.explained_variance_ratio_[1]*100, 2),figtype='pdf')

# get 2D biplot

cluster.biplot(cscore=pca_score, loadings=loadings, labels=loadings_df.index,
var1=round(pca.explained_variance_ratio_[0]*100, 2),
    var2=round(pca.explained_variance_ratio_[1]*100, 2))

# get correlation matrix plot for loadings

important_loding=loadings_df.iloc[:,0:6].copy()
important_loding=important_loding.round(2)
cut=1e-1;

tmpdf=abs(important_loding)>cut
important_loding=important_loding[tmpdf.any(axis='columns')]
# for col in important_loding:
#     important_loding[col] = important_loding[col].map('${:,.2f}'.format)


# ax = sns.clustermap(important_loding, annot=True, cmap='Spectral',col_cluster=False)
# plt.savefig('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/pca/Loading_pathways.pdf')

### get covariance matrix
pathCorr=np.corrcoef(pathway_fluxT.T)
# for i in range(pathCorr.shape[0]):
#     for j in range(pathCorr.shape[1]):
#         if i<j:
#             pathCorr[i,j]=0





df=pd.DataFrame(index=subsystems,columns=subsystems,data=pathCorr)

df=df.fillna(0)
mask = np.zeros_like(df.values)
mask[np.triu_indices_from(mask)] = True

with sns.axes_style("white"):
    f, ax = plt.subplots(figsize=(21, 21))
    ax = sns.heatmap(df,mask=mask,square=True,cmap='seismic')

# plt.savefig('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/pca/Pathway_pathway_correlation.pdf')
