getwd()
setwd("/Users/azamsaied/Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA")
#setwd("~/projects/MAGMA_Celltyping_benchmarking/")
#setwd("/home/rstudio/projects/MAGMA_Celltyping_benchmarking/")

#Lets install the MAGMA Celltypubg package
library(credentials)
Sys.getenv("GITHUB_PAT")
Sys.setenv(GITHUB_PAT = "xxxxxxxxxxx")
#set_github_pat()
#gitcreds::gitcreds_set(url = "https://github.com")
#gh::gh_token()


if(!require("remotes")) install.packages("remotes")
remotes::install_github("neurogenomics/MAGMA_Celltyping", force = TRUE)

#Lets install MungeSumStats
remotes::install_github("neurogenomics/MungeSumStats", force = TRUE)

#Install EWCE
BiocManager::install(version = '3.19')
BiocManager::install("rtracklayer")
BiocManager::install("piggyback")
BiocManager::install("EWCE", version = '3.19', force = TRUE) 

BiocManager::install(c("SNPlocs.Hsapiens.dbSNP144.GRCh37",
                       "SNPlocs.Hsapiens.dbSNP144.GRCh38",
                       "SNPlocs.Hsapiens.dbSNP155.GRCh37",
                       "SNPlocs.Hsapiens.dbSNP155.GRCh38",
                       "BSgenome.Hsapiens.1000genomes.hs37d5",
                       "BSgenome.Hsapiens.NCBI.GRCh38"))

library(MAGMA.Celltyping)
library(dplyr)
library(MungeSumstats)
library(EWCE)

#Setup - Specify where you want the large files to be downloaded to.
#NB - Set this folder to somewhere other than tempdir if you don't want the results to be deleted at the end.
#Since this is all just practice - I'm going to keep the tempdir
storage_dir <-("/Users/azamsaied/Documents/ResearchMagma Celltyping/Data")

##### Now lets try celltyping the munged sumstats ####

#Download GWAS summary statistics file to work on. This is pre-munged

# CTD = CellTypeDatasets 
# used as cell-type transcriptomic signature reference files
#Can make your own, or download an oven ready one
#ctd <- ewceData::ctd()
#ctd <- MAGMA.Celltyping::get_ctd("ctd_allKI")
des <- MAGMA.Celltyping::get_ctd("ctd_DescartesHuman")
#ctd <- get_ctd("ctd_Tasic")
#ctd <- get_ctd("ctd_DivSeq")
#ctd <- get_ctd("ctd_AIBS")
#dro <- get_ctd("ctd_DRONC_human")
#ctd <- get_ctd("ctd_DRONC_mouse")
#ctd <- get_ctd("ctd_BlueLake2018_FrontalCortexOnly")
#ctd <- get_ctd("ctd_BlueLake2018_VisualCortexOnly")
#ctd <- get_ctd("ctd_Saunders")


#Lets write some script forcing descartes ctd into different levels and saving it as such
des_2 <- des[c(2,3,4,5)]
names(des_2) <- c("level_1","level_2","level_3", "level_4")

des_3 <- des[c(3,4,5)]
names(des_3) <- c("level_1","level_2","level_3")
names(des_3) <- des[c(3,4,5)]

des_4 <- des[c(4,5)]
names(des_4) <- des[c(4,5)]

des_5 <- des[c(5)]
names(des_5) <- des[c(5)]

#Run cell-type enrichment analyses

#MAGMA.Celltyping has different functions for conducting various types of cell-type-specific enrichment tests on GWAS
#The celltype_associations_pipeline combines Linear enrichment, Top 10% enrichment, Conditional enrichment:

path_formatted <- "./Munged_cancer/colorectal_carcinoma_Auwerx_m.tsv.gz"
#path_formatted <- "./Munged_protein/CLPS_m.tsv.gz"

#Map SNPs to genes.
#Note - can input the genome build of the sum stats here - or it can be inferred if left NULL:
#Why not use 38?
genesOutPath <- MAGMA.Celltyping::map_snps_to_genes(path_formatted = path_formatted)
  
#Use thils one for Cancer
MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = "./Munged_cancer/MAGMA_Files/colorectal_carcinoma_Auwerx_m.tsv.35UP.10DOWN/",
  ctd = des_3,
  ctd_species = "human",   
  #ctd_levels = c(3),
  ctd_name = "ctd_Descartes_3", 
  run_linear = TRUE, 
  run_top10 = TRUE)

#Now lets plot the results 
#merge_results imports each of the MAGMA enrichment results files and merges them into one so that they can easily be plotted and further analysed.

merged_results <- MAGMA.Celltyping::merge_results(MAGMA_results = MAGMA_results)
knitr::kable(merged_results)

write.csv(merged_results, "./MAGMA_Celltype_Results/Descartes_L3/colorectal_carcinoma_Auwerx_L3_Merged_results.csv")
head(merged_results)
