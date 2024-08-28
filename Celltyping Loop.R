getwd()
setwd("/Users/azamsaied/Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA")
setwd("/Volumes/Expansion/R MAGMA/")
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
storage_dir <-("/Users/azamsaied/Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA/Data")

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
head(des_3)
names(des_2) <- c("level_1","level_2","level_3")

des_4 <- des[c(4,5)]
names(des_4) <- des[c(4,5)]

des_5 <- des[c(5)]
names(des_5) <- des[c(5)]

#Run cell-type enrichment analyses

#MAGMA.Celltyping has different functions for conducting various types of cell-type-specific enrichment tests on GWAS
#The celltype_associations_pipeline combines Linear enrichment, Top 10% enrichment, Conditional enrichment:

#setwd("/Users/azamsaied/Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA")
setwd("/Volumes/Expansion/R MAGMA/")
getwd()

#Re-write to continue loop if there is an error
pth <- "/Volumes/Expansion/R MAGMA/Munged_protein/m3000/"
ss <- list.files(pth)
#save_dir <- "/Users/azamsaied/Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA/CT_3000/"
save_dir <- "/Volumes/Expansion/R MAGMA/CT_3000/"

for(ss_i in ss) {
  print(ss_i)
  ss_i_ne <- strsplit(strsplit(ss_i, ".tsv")[[1]], "_")[[1]]
  prot_i <- ss_i_ne[[1]]
  syn_i <- ss_i_ne[[2]]
  magma_dir <- paste0("/Volumes/Expansion/R MAGMA/Munged_protein/m3000/MAGMA_Files/", prot_i, "_", syn_i, ".tsv.gz.35UP.10DOWN/")
  
  # Only run if not analysed before
  if(!file.exists (paste0("/Volumes/Expansion/R MAGMA/CT_3000/",prot_i, "_", syn_i, "_L2_Merged_results.csv"))){
    tryCatch({
      path_formatted <- paste0(pth, ss_i)
      genesOutPath <- MAGMA.Celltyping::map_snps_to_genes(path_formatted = path_formatted)
      MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
        magma_dirs = magma_dir,
        ctd = des_2,
        ctd_species = "human",   
        ctd_name = "ctd_Descartes_2", 
        run_linear = TRUE, 
        run_top10 = TRUE
      )
      merged_results <- MAGMA.Celltyping::merge_results(MAGMA_results = MAGMA_results)
      write.csv(merged_results, paste0("/Volumes/Expansion/R MAGMA/CT_3000/", prot_i, "_", syn_i, "_L2_Merged_results.csv"))
    }, error = function(e) {
      message("Error in processing: ", ss_i, " - ", e)
    })
  }
  
  # Delete sumstats
  #tryCatch({
    #file.remove(paste0(pth, ss_i))
  #}, error = function(e) {
    #message("Error in removing file: ", ss_i, " - ", e)
  #})
}





#Re-write to continue loop if there is an error
pth <- "/Volumes/Expansion/R MAGMA/Munged_protein/m3000/"
ss <- list.files(pth)
save_dir <- "/Users/azamsaied/Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA/CT_3000/"


for(ss_i in ss) {
  print(ss_i)
  ss_i_ne <- strsplit(strsplit(ss_i, ".tsv")[[1]], "_")[[1]]
  prot_i <- ss_i_ne[[1]]
  syn_i <- ss_i_ne[[2]]
  magma_dir <- paste0("/Volumes/Expansion/R MAGMA/Munged_protein/m3000/MAGMA_Files/", prot_i, "_", syn_i, ".tsv.gz.35UP.10DOWN/")
  
  # Only run if not analysed before
  if(!dir.exists(magma_dir)){
    tryCatch({
      path_formatted <- paste0(pth, ss_i)
      genesOutPath <- MAGMA.Celltyping::map_snps_to_genes(path_formatted = path_formatted)
      MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
        magma_dirs = magma_dir,
        ctd = des_3,
        ctd_species = "human",   
        ctd_name = "ctd_Descartes_3", 
        run_linear = TRUE, 
        run_top10 = TRUE
      )
      merged_results <- MAGMA.Celltyping::merge_results(MAGMA_results = MAGMA_results)
      write.csv(merged_results, paste0("/Volumes/Expansion/R MAGMA/CT_3000/", prot_i, "_", syn_i, "_L3_Merged_results.csv"))
    }, error = function(e) {
      message("Error in processing: ", ss_i, " - ", e)
    })
  }
  
  # Delete sumstats
  #tryCatch({
  #file.remove(paste0(pth, ss_i))
  #}, error = function(e) {
  #message("Error in removing file: ", ss_i, " - ", e)
  #})
}