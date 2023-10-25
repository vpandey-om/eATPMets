## Exploring the Impact of eATP on E.coli Metabolism Using REMI

We used [iML1515](http://bigg.ucsd.edu/models/iML1515)), a genome-scale metabolic model, to investigate the effects of eATP on E.coli. This investigation extended to determining the effect of eATP in two unique mediums, LB and M9. Our models for the M9 and LB mediums were carefully tweaked to match the particular compositions of these mediums, and these adjusted models can be found in the models folder.

We used [REMI](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007036) to effectively combine perturbation data linked to metabolites and genes. We employed dedicated scripts for constructing REMI models to acquire insights into the influence of eATP in the setting of both LB and M9 mediums:

1. [LB_medium_REMI_script.m](https://github.com/vpandey-om/eATPMets/blob/main/codes/matlabscript/LB_medium_REMI_script.m)
2. [M9_medium_REMI_script.m](https://github.com/vpandey-om/eATPMets/blob/main/codes/matlabscript/M9_medium_REMI_script.m)

