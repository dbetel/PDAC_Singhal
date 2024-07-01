# PDAC_Singhal
Code repository for Singhal et al. PDAC study

Analysis of scRNA-seq data in python


## Analysis Notebooks

### ScRNA-seq Data 

Each of the 4 datasets: 1 mouse and 3 PDXs (1,GA60, 106) have two notebooks for analysis a quality control (QC) notebook and analysis notebook.

### GSEA

Gene set enrichment analysis to determine programs upregulated in basal and classical cell states before and after treatment.Package details can be found here: [gseapy](https://gseapy.readthedocs.io/en/latest/introduction.html).


### inferCNV

Copy number variation inference to infer clonality of basal and classical cell states before and after treatment. Package deatils can be found here: [inferCNVpy](https://infercnvpy.readthedocs.io/en/latest/)

## Data

Data used for Singhal etal can be found at GEOXXXXX.

Gene position files for running inferCNV.

Conda environment to run notebooks can be created using singhalEtal.yml file. 

## Figures

Figures output folder