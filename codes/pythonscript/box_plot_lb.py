import pandas as pd

import plotnine

from plotnine import *

colorDict={'LB':"#F59751",'LB+ATP':"#7BB4BE"}


# df2=df[~(df['compound']=='TKT1')]

rxns=['GALUi','UMPK','INSH','TRPtipp']

for rxn in rxns:
    # rxn=''
    filename='/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/lb_rxns_flux_plot_'+rxn+'.txt'
    pdffile='/Users/vikash/Documents/MATLAB/eatp_metabolism/results/modelSim/lb_rxns_flux_plot_'+rxn+'.pdf'
    df=pd.read_csv(filename,sep='\t')
    manufacturer_cat = pd.Categorical(df['condi'], categories=['LB','LB+ATP'])
    df = df.assign(condition = manufacturer_cat)
    p10 = (
        ggplot(df, aes("compound", "wtFlux", fill="condition"))
        + geom_boxplot(outlier_shape = '')
        + xlab("Reaction")
        + geom_point(aes(fill = "condition"), size = 2, position = position_jitterdodge())
        # + facet_grid("compound ~ .", scales = "free_y")
        + ylab("Flux (mmol/l/h)")
        + scale_fill_manual(colorDict)
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
        )
        # + scale_fill_brewer(type="qual", palette="Accent",
        #                     name="Condition")
    )

    p10.save(pdffile)
