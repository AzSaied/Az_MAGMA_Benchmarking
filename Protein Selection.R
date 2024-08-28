setwd("/Volumes/Expansion/R MAGMA")
getwd()
# Load required libraries
library(dplyr)
library(readr)

### Step 1 - Prepping the data

# Read the data
p_atlas <- read_tsv("./rna_single_cell_type_tissue.tsv")

# Display the first few rows and a specific row
head(p_atlas)
p_atlas[1, ]

#We now have a dataframe with 11,185,674 rows and 7 columns
# Gene, Gene name, Tissue, Cluster, Cell type, Read count, nTPM
#The naming of the columns is a little confusing. I am going to rename 'Gene name' to 'Protein'
#I am aware these are not the same thing and there might not be a perfect correlation between the two - but the protein atlas data I have downloaded seems to measure RNA reads of the gene pre translation into protein. In the absense of proteomic measurements, it seems reasonable to use the read count for RNA transcripts of a gene, as a proxy for expression of a protein within a given cell type...


# Rename the 'Gene name' column to 'Protein'
p_atlas <- p_atlas %>% rename(Protein = `Gene name`)

# I find the 'gene' column unhelpful and annoying. I am going to get rid of it for now
# Remove the 'Gene' column and set 'Protein' as the index
p_atlas <- p_atlas %>% select(Protein, Tissue, Cluster, `Cell type`, `Read count`, nTPM)

# Display the first few rows
head(p_atlas)

#Lets see how many unique proteins we have

# Count unique proteins
unique_proteins <- p_atlas %>% distinct(Protein)
n_distinct(p_atlas$Protein)

#Wonderful - we now have 20,072 unique proteins

# Count occurrences of each protein
p_atlas %>% count(Protein)

#I notice that the counts for each goes in multiples of 557 (557, 1114, or 1671)
#I suspect (/hope) that this means there are 557 cell types in the protein atlas, though why some proteins appear two or three times per cell type isn't clear to me yet

#Now I want to remove any proteins that I do not have GWAS data for

#Let's look at the synapse UKB-PPP data
p_GWAS <- read_tsv("./olink_protein_map_3k_v1.tsv", show_col_types = FALSE)

#Here we have a list of 2,923 GWAS which we have good proteomic data for

# Display the first few rows and a specific row
head(p_GWAS)
p_GWAS[1, ]

# Rename the 'Assay' column to 'Protein'
p_GWAS <- p_GWAS %>% rename(Protein = `Assay`)

# Remove duplicate proteins
p_GWAS <- p_GWAS %>% distinct(Protein, .keep_all = TRUE)

#So 2,923 of the UKB:PPP gwas list had unique proteins

#Lastly I want to remove any proteins that are in the protein atlas but that we do not have GWAS data for:
p_atlas_final <- p_atlas %>% filter(Protein %in% p_GWAS$Protein)

#Wonderful!
#We now have a lovely df which has 1,609,173 cell measurements and 8 columns

#I think i'd find the next step easier if my column headings were 1 word
p_atlas_final <- p_atlas_final %>% rename(Cell_type = `Cell type`, Read_count = `Read count`)
head(p_atlas_final)


### Step 2 - Data exploration
#I want to explore the data - look at proteins that are expressed in one cell type but not a (closely related) cell type.


# Filter data for a specific protein
p_atlas_final %>% filter(Protein == 'FXN')

# Count occurrences of each tissue
p_atlas_final %>% count(Tissue)

#So there are 31 Tissues.
#brain, fallopian tube, placenta, lung, breast, endometrium, salivary gland, adipose tissue, tongue, testis, 
#vascular, thymus, skeletal muscle, liver, skin, pancreas, ovary, esophagus, prostate, bronchus, kidney, 
#small intestine, eye, pbmc, colon, bone marrow, rectum, spleen, stomach, heart muscle, lymph node)



#This is a problem. I am looking for discrepancy BETWEEN different cell types, but if there is discrepance WITHIN a given cell type - how can I interpret the discrepancy BETWEEN cells?

#One option is to take the MAX value. If the protein can be expressed in a cell, it can be expressed in a cell. End of.
#If one cardiomyocyte does not express a protein, but another does - the fact remains that the protein is being expressed in cardiomyocytes.
#I am starting with 557 cell types
#Lets calculate how many there will be if I remove all duplicate cell types:

# Count unique cell types
unique_cell_types <- p_atlas_final %>% distinct(Cell_type)
n_distinct(p_atlas_final$Cell_type)

# Count occurrences of each cell type
p_atlas_final %>% count(Cell_type)

#This takes me down to 85 cell types.
#This shows that 25 cell types are unique, 10 appear twice, 12 appear three times, 6 appear four times, 32 appear five or more times

#Lets remove duplicates keeping those with the highest value.
#I want to be smart about this, I only want a given cell type to appear once for each tissue type
#I want the read count to reflect the highest value for a given cell type - regardless of whether that count came from another tissue.
#eg if a T-cell in the LN expressed a high count, even when looking at T-cells in the bone marrow, I know the cells are capable of expressing that RNA (and hence protein)

#This is a slightly trickier way of doing it - but let's try
#I want to look at every occourance of a cell type, across all tissues and find the max value for the read count. 
#I then want to apply that read count to every one of those cell types across all tissues.
#I then want to remove duplicates - BUT ONLY WITHIN A GIVEN TISSUE
#I think groupby will be useful here


# Group by 'Protein' and 'Cell_type', then find max 'Read_count'
tissue_ignore_df <- p_atlas_final %>% group_by(Protein, Cell_type) %>% summarize(Max_Read_count = max(Read_count)) %>% ungroup()

# Merge with original data
p_atlas_final <- p_atlas_final %>% left_join(tissue_ignore_df, by = c("Protein", "Cell_type"))

# Calculate difference in read counts
p_atlas_final <- p_atlas_final %>% mutate(RC_diff = Read_count - Max_Read_count)

#Ok great - I now have a tidy dataframe which expresses the read count for a given cell ina given tissue, but also the max read count of that cell anywhere in the body
#This is useful in determining whether that cell is capable of expressing a given protein

#p_atlas_final_V2 <- p_atlas_final

# Now lets write some code which lets us look at a specific tissue type - to help us choose proteins which distinguish between cell types within a given tissue








###That's great!###


#Now lets do it for another tissue - Lung
lung_atlas <- p_atlas_final %>% filter(Tissue == 'lung')

#Because I am only interested in proteins which can be used to form a 'True Negative' test, I need at least one cell type within the tissue in question, to have zero max expression of that protein
#The following lines of code will filter out proteins which have no zero expression

# Filter skin_atlas to keep only proteins with some 'zero' expression
lung_0_protein <- lung_atlas %>% filter(Max_Read_count == 0) %>% distinct(Protein)
lung_atlas_zerofilter <- lung_atlas %>% filter(Protein %in% lung_0_protein$Protein)

#We can then write it to a CSV to peruse it in excel, or just look at it as an object in R-studio
write_csv(lung_atlas_zerofilter, 'lung_atlas_zerofilter.csv')

#Lets try being more selective about which cells types have a max read count of 0
dummy <- BM_atlas %>% filter(Cell_type != 'erythroid cells' & Cell_type != 'plasma cells')
BM_0_protein <- dummy %>% filter(Max_Read_count == 0) %>% distinct(Protein)
BM_atlas_zerofilter <- dummy %>% filter(Protein %in% BM_0_protein$Protein)
write_csv(BM_atlas_zerofilter, 'BM_atlas_zerofilter_v2.csv')




#OK - Now I want to just create a list of all of the proteins which will eventially be worth using

#Lets try an pull off a list of proteins which have at least 1 zero count anywhere:
zero_proteins <- p_atlas_final %>% filter(Max_Read_count == 0) %>% distinct(Protein)
total_p_atlas_0 <- p_atlas_final %>% filter(Protein %in% zero_proteins$Protein)

#Lets filter out anything which isn't a zero count
zero_proteins_2 <- p_atlas_final %>% filter(Max_Read_count == 0)
#And then make a table counting how many 'zeros' each protein has
counts <- as.data.frame(table(zero_proteins_2['Protein']))
write_csv(counts, 'PPP_counts.csv')

dummy <- total_p_atlas_0 %>% filter(Protein == 'ALPI')

#I want to create a column which concat's the tissue and cell type columns
library(tidyr)
total_p_atlas_0 = total_p_atlas_0 %>% 
  unite(combined, Tissue, Cell_type, sep = "_ ", remove = FALSE)

#Now I want to read in my list of final tissue_celltypes for the truth matrix
TM.df <- read_csv("TM_Tissue_Celltypes.csv")

#Lets filter out all the the celltypes we are not keeping in our final TM
p_atlas_final_celltypes <- merge(x = total_p_atlas_0, y = TM.df, by = c("combined"))
p_atlas_final_celltypes <- p_atlas_final_celltypes %>% select(-Tissue, -Cluster, -Cell_type, -Read_count, -nTPM, -RC_diff)

#Lets save this a a CSV file...but then start working on it
write_csv(p_atlas_final_celltypes, "/Volumes/Expansion/R MAGMA/HPA_Celltyes_Proteins_Counts.csv")



##NEGATIVE TRUE

#Lets work on making the datatable I want- starting with the p_atlas_final_celltypes
#I want to translate the cell types from HPA into Descartes
#Lets read in the key CSV file
HPA_Des3_key <- read_csv("./Keys/HPA-Descartes_L3.csv")

#Now lets try to read and replace from the HPA Des3 Key
library(dplyr)
p_atlas_final_celltypes$Des3 <- ifelse(p_atlas_final_celltypes$combined %in% HPA_Des3_key$HPA, 
                          HPA_Des3_key$`Descartes L3`[match(p_atlas_final_celltypes$combined, HPA_Des3_key$HPA)], 
                          p_atlas_final_celltypes$combined[match(p_atlas_final_celltypes$combined, HPA_Des3_key$`Descartes L3`)])
p_atlas_final_celltypes <- p_atlas_final_celltypes[!duplicated(p_atlas_final_celltypes), ]

#Great - now we have a lovely table with the translated cell types for descartes L3.
#Now lets filter out all of the non 0 max values
Noughty_list <- p_atlas_final_celltypes %>% filter(Max_Read_count == '0')
#Lets just keep the proteins in my analysis
#Noughty_list <- p_atlas_final_celltypes %>% filter(Protein == c['FSHB','CABP2','FGF6',	'PTH',	'CEACAM18',	'IL3',	'OMP',	'IFNW1',	'MAGEA3',	'ALPI'])
#CT_protein_list <- read.csv(file=#path to your csv within quote marks "myfile.csv")
#For now lets just use: Positive_proteins
library(tidyverse)
Noughty_list <- filter(Noughty_list,Protein %in% Positive_proteins$Var1)

#I need this to turn into an array with 2 columns I assume
#lets concatinate the protein and the Des3 cell type name
Noughty_list$P_CT <- paste0(Noughty_list$Protein,"_",Noughty_list$Des3)

#Now lets just grab the columns I need to make the array that PRAUC needs
PRAUC_Negative_True <- Noughty_list %>% select(P_CT, Max_Read_count)















#NEGATIVE PREDICTION
#(Using the poritive porotein list)

# Load required libraries
library(dplyr)
library(readr)
library(tidyverse)

# Load protein and associated synapse keys
Synapse_key <- data.table::fread("./Synapse_keys.csv")
Synapse_key[, prot := sub("\\_.*", "", Name)]

# Load the list of proteins of interest
Positive_proteins <- read_csv("/Volumes/Expansion/R MAGMA/Pos_prot_2.csv")

# Initialize the PM dataframe with the list of cell types
PM <- read_csv("/Volumes/Expansion/R MAGMA/Protein_matrix.csv")

# Define the path and get list of files
pth <- "/Volumes/Expansion/R MAGMA/CT_3000/"
ss <- list.files(pth, pattern = "*.csv")
print(ss)
# Debugging print
print("Initial PM structure:")
print(head(PM))

# Loop through each file in the directory
for(ss_i in ss) {
  print(paste("Processing file:", ss_i))
  ss_i_ne <- strsplit(strsplit(ss_i, ".csv")[[1]], "_")[[1]]
  prot_i <- ss_i_ne[[1]]
  syn_i <- ss_i_ne[[2]]
  
  if(prot_i %in% Positive_proteins$Var1) {
    tryCatch({
      # Read the CSV file
      print(prot_i)
      temp <- read_csv(paste0(pth, ss_i))
      temp <- temp %>% select(Celltype, FDR)
      temp <- temp %>% 
        group_by(Celltype) %>% 
        summarise(FDR = min(FDR))
      
      # Rename the FDR column to the current protein name
      temp <- temp %>% rename(!!prot_i := 'FDR')
      
      # Debugging print
      print(paste("Merging data for protein:", prot_i))
      print(head(temp))
      
      # Merge PM and temp by Celltype
      PM <- PM %>% left_join(temp, by = "Celltype")
      
      # Debugging print to check the structure after merge
      print("PM structure after merge:")
      print(head(PM))
      
    }, error = function(e) {
      print(paste0("Error in ", prot_i))
    })
  }
}

# Final structure check
print("Final PM structure:")
print(head(PM))

# Save the final PM dataframe to a CSV file if needed
write_csv(PM, "/Volumes/Expansion/R MAGMA/Updated_Protein_matrix.csv")



































#Now lets make mini dataframes lookint at one protein at a time
FSHB <- p_atlas_final_celltypes %>% filter(Protein == 'FSHB')
FSHB <- FSHB[!duplicated(FSHB), ]
FSHB <- FSHB %>% select(-Protein)
FSHB <- FSHB %>% rename(FSHB = `Max_Read_count`)
CABP2 <- p_atlas_final_celltypes %>% filter(Protein == 'CABP2')
CABP2 <- CABP2[!duplicated(CABP2), ]
CABP2 <- CABP2 %>% select(-Protein)
CABP2 <- CABP2 %>% rename(CABP2 = `Max_Read_count`)
FGF6 <- p_atlas_final_celltypes %>% filter(Protein == 'FGF6')
FGF6 <- FGF6[!duplicated(FGF6), ]
FGF6 <- FGF6 %>% select(-Protein)
FGF6 <- FGF6 %>% rename(FGF6 = `Max_Read_count`)
PTH <- p_atlas_final_celltypes %>% filter(Protein == 'PTH')
PTH <- PTH[!duplicated(PTH), ]
PTH <- PTH %>% select(-Protein)
PTH <- PTH %>% rename(PTH = `Max_Read_count`)
MLN <- p_atlas_final_celltypes %>% filter(Protein == 'MLN')
MLN <- MLN[!duplicated(MLN), ]
MLN <- MLN %>% select(-Protein)
MLN <- MLN %>% rename(MLN = `Max_Read_count`)
CEACAM18 <- p_atlas_final_celltypes %>% filter(Protein == 'CEACAM18')
CEACAM18 <- CEACAM18[!duplicated(CEACAM18), ]
CEACAM18 <- CEACAM18 %>% select(-Protein)
CEACAM18 <- CEACAM18 %>% rename(CEACAM18 = `Max_Read_count`)
IL3 <- p_atlas_final_celltypes %>% filter(Protein == 'IL3')
IL3 <- IL3[!duplicated(IL3), ]
IL3 <- IL3 %>% select(-Protein)
IL3 <- IL3 %>% rename(IL3 = `Max_Read_count`)
OMP <- p_atlas_final_celltypes %>% filter(Protein == 'OMP')
OMP <- OMP[!duplicated(OMP), ]
OMP <- OMP %>% select(-Protein)
OMP <- OMP %>% rename(OMP = `Max_Read_count`)
IFNW1 <- p_atlas_final_celltypes %>% filter(Protein == 'IFNW1')
IFNW1 <- IFNW1[!duplicated(IFNW1), ]
IFNW1 <- IFNW1 %>% select(-Protein)
IFNW1 <- IFNW1 %>% rename(IFNW1 = `Max_Read_count`)
MAGEA3 <- p_atlas_final_celltypes %>% filter(Protein == 'MAGEA3')
MAGEA3 <- MAGEA3[!duplicated(MAGEA3), ]
MAGEA3 <- MAGEA3 %>% select(-Protein)
MAGEA3 <- MAGEA3 %>% rename(MAGEA3 = `Max_Read_count`)

#Lets add ALPI to the list
ALPI <- p_atlas_final_celltypes %>% filter(Protein == 'ALPI')
ALPI <- ALPI[!duplicated(ALPI), ]
ALPI <- ALPI %>% select(-Protein)
ALPI <- ALPI %>% rename(ALPI = `Max_Read_count`)

#Now lets combine them into the Truth Matrix
#TM <- merge(TM.df, FSHB, CABP2, FGF6, PTH, MLN, CEACAM18, IL3, OMP, IFNW1, MAGEA3)
TM <- merge(TM.df, FSHB)
TM <- merge(TM, CABP2)
TM <- merge(TM, FGF6)
TM <- merge(TM, PTH)
TM <- merge(TM, MLN)
TM <- merge(TM, CEACAM18)
TM <- merge(TM, IL3)
TM <- merge(TM, OMP)
TM <- merge(TM, IFNW1)
TM <- merge(TM, MAGEA3)
TM <- merge(TM, ALPI)

#LEts export that into an excel file to have a play with it
write_csv(TM, 'TM_10p.csv')


#Lets look at some of our shortlisted proteins
dummy <- total_p_atlas_0 %>% filter(Protein == 'FGF6')
