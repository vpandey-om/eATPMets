import pandas as pd
import sys,os
from scipy.stats import gmean
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt

# setting path
# getting the name of the directory
# where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))
# Getting the parent directory name
# where the current directory is present.
eatp_metabolism = os.path.dirname(current)
# adding the parent directory to
# the sys.path.

df_1m60=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data','Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_1mMATP_60min_vikash.xlsx'))
df_1m180=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data','Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_1mMATP_180min_vikash.xlsx'))
df_100u60=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data','Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_100uM_ATP_60min_Vikash.xlsx'))
df_500u60=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data','Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_500uM_ATP_60min_Vikash.xlsx'))
df_500u180=pd.read_excel(os.path.join(eatp_metabolism,'Metabolomic_data','Bacterial_pellet_in_Rich_medium_(LB)_-_FC_analysis_', 'LB_500uM_ATP_180min_Vikash.xlsx'))

dfList=[df_1m60,df_1m180,df_100u60,df_500u60,df_500u180]
res_df=df_1m60.merge(df_1m180, on='Column1',suffixes=['_1m60','_1m180'],how='outer')
res_df=res_df.merge(df_100u60, on='Column1',suffixes=['','_100u60'],how='outer')
res_df=res_df.merge(df_500u60, on='Column1',suffixes=['','_500u60'],how='outer')
res_df=res_df.merge(df_500u180, on='Column1',suffixes=['','_500u180'],how='outer')


def plot_scatter(ax,df,title_text):
    x=df['log2(FC)'].values
    y= df['-log10(p)'].values
    xmax=np.max(x)
    ymax=np.max(x)
    ax.scatter(x,y)
    for i, txt in enumerate(df['Column1'].values):
        ax.annotate(txt, (x[i]-0.05, y[i]-0.05),fontsize = 6)
    ax.set_xlabel('log2(FC)')
    ax.set_ylabel('-log10(p)')
    ax.set_title(title_text)
    if xmax<2.5:
        ax.set_xlim(-2.5,2.5)
    if xmax<2:
        ax.set_ylim(0,2)



## plot scatter
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20*0.8, 30*0.8))
plot_scatter(axes[0][1],df_100u60,'100uM60min')
plot_scatter(axes[1][1],df_500u60,'500uM60min')
plot_scatter(axes[1][0],df_500u180,'500uM180min')
plot_scatter(axes[2][1],df_1m60,'1mM60min')
plot_scatter(axes[2][0],df_1m180,'1mM180min')


fig.savefig('all_rich_data.pdf')

# axes[0][0].scatter(df_1m60['log2(FC)'].values, df_1m60['-log10(p)'].values)
# axes[0][0].set_xlabel('log2(FC)')
# axes[0][0].set_ylabel('-log10(p)')

# axes[0][1].scatter(df_1m180['log2(FC)'].values, df_1m180['-log10(p)'].values)
# axes[1][0].scatter(df_100u60['log2(FC)'].values, df_100u60['-log10(p)'].values)
# axes[1][1].scatter(df_500u60['log2(FC)'].values, df_500u60['-log10(p)'].values)
# axes[2][0].scatter(df_500u180['log2(FC)'].values, df_500u180['-log10(p)'].values)





import pdb; pdb.set_trace()
