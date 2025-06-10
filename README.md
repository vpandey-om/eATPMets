## Exploring the Impact of eATP on _E. coli_ Metabolism Using REMI

We used [iML1515](http://bigg.ucsd.edu/models/iML1515), a genome-scale metabolic model, to investigate the effects of eATP on _E. coli_. We investigated the effect of eATP in two unique media, rich medium (LB) and minimal medium (M9). Our models for the M9 and LB media were carefully tweaked to match the particular compositions of these media, and these adjusted models can be found in the models folder.

We used [REMI](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007036) to effectively combine perturbation data linked to metabolites and genes. We employed dedicated scripts for constructing REMI models to acquire insights into the influence of eATP in the setting of both LB and M9 mediums:

1. [LB_medium_REMI_script.m](https://github.com/vpandey-om/eATPMets/blob/main/codes/matlabscript/LB_medium_REMI_script.m)
2. [M9_medium_REMI_script.m](https://github.com/vpandey-om/eATPMets/blob/main/codes/matlabscript/M9_medium_REMI_script.m)


# NCA & Clustering Pipeline (`eATPMets/nca`)

This guide outlines a streamlined workflow for performing Network Component Analysis (NCA) using CRP-regulated genes and clustering the results with Scanpy.

---

## 🔹 Steps Overview

1. **Collect CRP-Regulated Genes**  
   Download CRP target genes from [RegulonDB](https://regulondb.ccg.unam.mx/) and save them in CSV or TSV format.

2. **Preprocess Data**  
   Run the preprocessing Python script:
   ```bash
   python scratch_fastnca.py

### 3. Run NCA and Cluster the Results

First, perform NCA using the MATLAB scripts:




```matlab
performNCA
fastNCA

Then, cluster tusing the following Python script:
```python python_scripts/clustering.py


