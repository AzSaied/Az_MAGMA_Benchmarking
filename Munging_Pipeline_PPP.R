getwd()
setwd("/home/rstudio/projects/as6820/")
#default is "/home/rstudio/projects/MAGMA_Celltyping_benchmarking"
options(timeout = max(10000, getOption("timeout")))
#Let's install all the packages we'll need to start:

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("MungeSumstats")

#BiocManager::install(c("SNPlocs.Hsapiens.dbSNP144.GRCh37",
                       "SNPlocs.Hsapiens.dbSNP144.GRCh38",
                       "SNPlocs.Hsapiens.dbSNP155.GRCh37",
                       "SNPlocs.Hsapiens.dbSNP155.GRCh38",
                       "BSgenome.Hsapiens.1000genomes.hs37d5",
                       "BSgenome.Hsapiens.NCBI.GRCh38"))

#if(!require("remotes")) install.packages("remotes")
#remotes::install_github("neurogenomics/MAGMA_Celltyping")

#install.packages("data.table")
#install.packages("gitcreds")
#library(gitcreds)
#gitcreds_set()

#Load our librarys
library(MAGMA.Celltyping)
library(data.table)
library(dplyr)
library(MungeSumstats)
library(R.utils)

#I want the files to be stored in my working directory - so I will keep this hashed out
#storage_dir <- tempdir("/home/rstudio/projects/as6820")


#I Want all of this to be able to run using MungeSS's abilty to search for datasets
metagwas <- MungeSumstats::find_sumstats(authors = "Walters RG")
head(metagwas,100)
ids <- (dplyr::arrange(metagwas, nsnp))$id  

datasets <- MungeSumstats::import_sumstats(ids = "ebi-a-GCST90013410",
                                           ref_genome = "GRCH37")

#Use this line of code to download a tsv file and safe it locally for munging
download.file('http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018838/GCST90018838_buildGRCh37.tsv.gz', 
              destfile = "Endometrial_Sakaue.tsv.gz", 
              method = "wget", extra = "-r -p --random-wait")

#Or if using Synapse - use this:
#install.packages("reticulate")
#install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapser)
#install.packages("synapserutils", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapserutils) 
synLogin(authToken="XXXXXXXX") 

# Obtain a pointer and download the data 
#syn52361407 < - synGet(entity=' syn52361407 ') 

#Lets try downloading and munging protein FSHB
# Ensemble ID ENSG00000170465
GAGE2A <- synGet(entity ='syn52362803')

#clicking on the FSHB R object this creates, provides a path to where the TSV file is saved:
#  "/home/rstudio/.synapseCache/580/131009580/"

#I'm going to try to untar the files:
help("archive_extract")
archive_extract("/home/rstudio/.synapseCache/580/131009580/GAGE2A_Q6NT46_OID31199_v1_Oncology_II.tar", "/home/rstudio/projects/as6820/GAGE2A")
#untar("/home/rstudio/projects/as6820/GAGE2A/GAGE2A_Q6NT46_OID31199_v1_Oncology_II/discovery_chr1_GAGE2A:Q6NT46:OID31199:v1:Oncology_II.gz")
#Now to combine all of the separate TAR files into one big one

setwd("/home/rstudio/projects/as6820/GAGE2A/GAGE2A_Q6NT46_OID31199_v1_Oncology_II/")
chrs <- paste0("chr",1:22)
all_chrs_dat <- vector(mode="list",length = length(chrs))
names(all_chrs_dat) <- chrs
for(chr_i in chrs){
  pth_i <- paste0("discovery_",chr_i,
                  "_GAGE2A:Q6NT46:OID31199:v1:Oncology_II.gz")
  tmp <- data.table::fread(pth_i)
  all_chrs_dat[[chr_i]] <- tmp
}

dat <- data.table::rbindlist(all_chrs_dat)
#save combined

setwd("/home/rstudio/projects/as6820/")

#Lets try addint my own column with a P value in it:
dat$P <- 0
dat$P <- 10^ (-1 * dat$LOG10P)
#and I need to help it recognise the SNP column. Lets remove the ':imp:v1' from the end of every cell in ID
dat$SNP = substr(dat$ID,1,nchar(dat$ID)-7)
dat2 <- select(dat, SNP, CHROM, GENPOS, ALLELE0, ALLELE1, TEST, INFO, N, BETA, SE, CHISQ, P)
names(dat2) <- c("Name", "chromosome", "base_pair_location", "other_allele", "effect_allele", "Model", "INFO", "N", "BETA", "SE", "CHISQ", "p_value")



##### Ok - let's try munging #####

##Set the path to the downloaded file##
setwd("/home/rstudio/projects/as6820/")
reformatted <- MungeSumstats::format_sumstats(path=dat2,
                                              ref_genome = "GRCh38",
                                              save_path = "./GAGE2A_m.tsv.gz")

#And that's the Munging complete

head(reformatted, 10)
