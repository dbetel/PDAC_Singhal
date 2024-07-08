# PDAC_Singhal
Code repository for Singhal et al. PDAC study


[A classical epithelial state drives acute resistance to KRAS inhibition in pancreas cancer](https://aacrjournals.org/cancerdiscovery/article/doi/10.1158/2159-8290.CD-24-0740/746288/A-classical-epithelial-state-drives-acute)

Analysis of scRNA-seq data in python


## Analysis Notebooks

### ScRNA-seq Data 

Each of the 4 datasets: 1 mouse and 3 PDXs (1,G A60, 106) have two notebooks for analysis a quality control (QC) notebook and analysis notebook.

### GSEA

Gene set enrichment analysis to determine programs upregulated in basal and classical cell states before and after treatment. Package details can be found here: [gseapy](https://gseapy.readthedocs.io/en/latest/introduction.html).


## Data

Data used for Singhal etal can be found at [GSE271300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271300).

Conda environment to run notebooks can be created using singhalEtal.yml file. 

## Figures

Figures output folder