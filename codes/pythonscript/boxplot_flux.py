import pandas as pd

import plotnine

from plotnine import *

colorDict={'M9':"#F59751",'M9+ATP':"#7BB4BE"}


# df2=df[~(df['compound']=='TKT1')]

rxns=['ACGK','ACGS','AMPN','CYTDK2','DADA','NNAM','PUNP2','PUNP6','TKT1','GLUPRT','GTHS']

for rxn in rxns:
    filename='/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/m9_rxns_flux_plot_'+rxn+'.txt'
    pdffile='/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/m9_rxns_flux_plot_'+rxn+'.pdf'
    df=pd.read_csv(filename,sep='\t')
    manufacturer_cat = pd.Categorical(df['condi'], categories=['M9','M9+ATP'])
    df = df.assign(condition = manufacturer_cat)
    p10 = (
        ggplot(df, aes("compound", "wtFlux", fill="condition"))
        + geom_boxplot(outlier_shape = '')
        + xlab("Reaction")
        + geom_point(aes(fill = "condition"), size = 2, position = position_jitterdodge())
        + scale_fill_manual(colorDict)
        # + scale_color_manual(colorDict)
        # + facet_grid("compound ~ .", scales = "free_y")
        + ylab("Flux (mmol/l/h)")
        + theme(
            legend_position="bottom",
            legend_direction="horizontal",
            legend_title_align="center",
            legend_box_spacing=0.6,
            legend_key=element_blank(),
            axis_line=element_line(size=1, colour="black"),
            panel_grid_major=element_line(colour="#d3d3d3"),
            panel_grid_minor=element_blank(),
            panel_border=element_blank(),
            panel_background=element_blank(),
            plot_title=element_text(size=15, family="Tahoma",
                                    face="bold"),
            text=element_text(family="Tahoma", size=11),
            axis_text_x=element_text(angle = 90,colour="black", size=10),
            axis_text_y=element_text(colour="black", size=10),
            # figure_size=(6.4, 12),
        )
        # + scale_fill_brewer(type="qual", palette="Accent",
        #                     name="Condition")
    )

    p10.save(pdffile)
