import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans,AffinityPropagation
import numpy as np
def draw_umap(n_neighbors=15, min_dist=0.1, n_components=3, metric='euclidean', title=''):
    fit = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric
    )
    u = fit.fit_transform(data);
    fig = plt.figure()
    if n_components == 1:
        ax = fig.add_subplot(111)
        ax.scatter(u[:,0], range(len(u)), c=colors,s=10)
    if n_components == 2:
        ax = fig.add_subplot(111)
        ax.scatter(u[:,0], u[:,1], c=colors,s=10)
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(u[:,0], u[:,1], u[:,2], c=colors, s=10)
    plt.title(title, fontsize=18)
    plt.savefig('testing'+str(n_neighbors)+'.pdf')

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


pca = PCA(n_components=2, svd_solver='full')

df=pd.read_csv('/Users/vikash/Documents/MATLAB/eatp_metabolism/results/LB_REMI.txt',header=None,sep='\t')
data=df.T.values
### get diffrence in flux 
y_pred = KMeans(n_clusters=20, random_state=0).fit_predict(data)

eatp_cluster=[]
ctl_cluster=[]
for c in range(20):
    xx=np.where(y_pred ==c)
    # print (len(xx[0]),xx[0])
    A_1=xx[0]
    cut1=(A_1 > 501).sum()/len(xx[0])
    cut2=(A_1 < 501).sum()/len(xx[0])

    if cut1>0.9:
        eatp_cluster.append(c)
    if cut2>0.9:
        ctl_cluster.append(c)


    print('cluster',c+1)
    print (len(xx[0]),((A_1 > 501).sum() == A_1.size).astype(np.int),(A_1 < 501).sum()/len(xx[0]),((A_1 < 501).sum() == A_1.size).astype(np.int))



eatpidx=[]
for c in eatp_cluster:
        xx=np.where(y_pred ==c)
        eatpidx.append(xx[0])
eatpidx = np.concatenate(eatpidx).ravel()
ctlIdx=[]
for c in ctl_cluster:
    xx=np.where(y_pred ==c)
    ctlIdx.append(xx[0])

ctlIdx = np.concatenate(ctlIdx).ravel()

print((eatpidx > 501).sum()/len(eatpidx))
print((ctlIdx  <501).sum()/len(ctlIdx ))
idx = np.concatenate([ctlIdx,eatpidx]).ravel()

data_filter=data[idx,:]

# af = AffinityPropagation(preference=-50, random_state=0).fit(data)
# cluster_centers_indices = af.cluster_centers_indices_
# labels = af.labels_

# n_clusters_ = len(cluster_centers_indices)



dist=pairwise_distances(data)
colors=[]
idx=[]
for i in range(len(ctlIdx)):
    colors.append('green')
    idx.append(i+1)
for i in range(len(eatpidx)):
    colors.append('red')
    idx.append(i+1)

pca.fit(data)
X_pca = pca.transform(data_filter)
draw_pca(X_pca,n_components=2)
import pdb; pdb.set_trace()




# for n in (2, 5, 10, 20, 50, 100, 200):
#     draw_umap(n_neighbors=n, title='n_neighbors = {}'.format(n))
