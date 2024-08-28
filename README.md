This is the readme file for my benchmarking of MAGMA_Celltyping project.

We downloaded accessed UKB-PPP GWAS data from https://www.synapse.org. The R-code can search and access the database directly, but needs the user to have an auth token from
the synapse website.

We downloaded HPA data from https://www.proteinatlas.org. 
We searched and downloaded non-proteomic GWAS summary statistics files from the opensource GWAS catalogue, and IEU open GWAS project. 

R version 4.4.0 and BiocManager version 1.30.23 were used to analyse this data to build the truth matrices used in our benchmarking. 

We used the opensource R package ‘MungeSumstats’  version 1.13.2 to convert the wide variety of downloaded GWAS summary statistics files 
into a standard format that could be run by the MAGMA_Celltyping package. 

MAGMA_Celltyping version 2.0.12, was used to perform the cell-typing using the Descartes L2 and L3 CTDs
