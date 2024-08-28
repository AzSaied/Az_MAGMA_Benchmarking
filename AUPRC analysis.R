setwd("/Volumes/Expansion/R MAGMA")
setwd("Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA/")
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

#OK - Now I want to just create a list of all of the proteins which will eventially be worth using

#Lets try an pull off a list of proteins which have at least 1 zero count anywhere:
zero_proteins <- p_atlas_final %>% filter(Max_Read_count == 0) %>% distinct(Protein)
total_p_atlas_0 <- p_atlas_final %>% filter(Protein %in% zero_proteins$Protein)

#Lets filter out anything which isn't a zero count
zero_proteins_2 <- p_atlas_final %>% filter(Max_Read_count == 0)
#And then make a table counting how many 'zeros' each protein has
counts <- as.data.frame(table(zero_proteins_2['Protein']))
#write_csv(counts, 'PPP_counts.csv')

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
#write_csv(p_atlas_final_celltypes, "./HPA_Celltyes_Proteins_Counts.csv")

##NEGATIVE TRUE
p_atlas_final_celltypes <- read_csv("./HPA_Celltyes_Proteins_Counts.csv")
head(p_atlas_final_celltypes)


### DESCARTES L3 ###

#Lets work on making the datatable I want- starting with the p_atlas_final_celltypes
#I want to translate the cell types from HPA into Descartes
#Lets read in the key CSV file
HPA_Des3_key <- read_csv("./Keys/HPA-Descartes_L3.csv")
HPA_Des2_key <- read_csv("./Keys/HPA-Descartes_L2.csv")
colnames(HPA_Des2_key) <- c("HPA", "Descartes L2")

#Now lets try to read and replace from the HPA Des3 Key
library(dplyr)
p_atlas_final_celltypes$Des3 <- ifelse(p_atlas_final_celltypes$combined %in% HPA_Des3_key$HPA, 
                                       HPA_Des3_key$`Descartes L3`[match(p_atlas_final_celltypes$combined, HPA_Des3_key$HPA)], 
                                       p_atlas_final_celltypes$combined[match(p_atlas_final_celltypes$combined, HPA_Des3_key$`Descartes L3`)])
p_atlas_final_celltypes <- p_atlas_final_celltypes[!duplicated(p_atlas_final_celltypes), ]
head(p_atlas_final_celltypes)

#Lets try adding a column with Des2
p_atlas_final_celltypes$Des2 <- ifelse(p_atlas_final_celltypes$combined %in% HPA_Des2_key$HPA, 
                                       HPA_Des2_key$`Descartes L2`[match(p_atlas_final_celltypes$combined, HPA_Des2_key$HPA)], 
                                       p_atlas_final_celltypes$combined[match(p_atlas_final_celltypes$combined, HPA_Des2_key$`Descartes L2`)])
p_atlas_final_celltypes <- p_atlas_final_celltypes[!duplicated(p_atlas_final_celltypes), ]
head(p_atlas_final_celltypes)

#Great - now we have a lovely table with the translated cell types for descartes L3 AND L2.
#Now lets filter out all of the non 0 max values
Noughty_list <- p_atlas_final_celltypes %>% filter(Max_Read_count == '0')
#Lets just keep the proteins in my analysis
#Noughty_list <- p_atlas_final_celltypes %>% filter(Protein == c['FSHB','CABP2','FGF6',	'PTH',	'CEACAM18',	'IL3',	'OMP',	'IFNW1',	'MAGEA3',	'ALPI'])
#CT_protein_list <- read.csv(file=#path to your csv within quote marks "myfile.csv")
#For now lets just use: Positive_proteins
library(tidyverse)
Positive_proteins <- read.csv("Pos_prot_2.csv")
Noughty_list <- filter(Noughty_list,Protein %in% Positive_proteins$Var1)
write_csv(Noughty_list, './Naughty_list.csv')


#I need this to turn into an array with 2 columns I assume
#lets concatinate the protein and the Des3 cell type name
Noughty_list$P_CT_L3 <- paste0(Noughty_list$Protein,"_",Noughty_list$Des3)
Noughty_list$P_CT_L2 <- paste0(Noughty_list$Protein,"_",Noughty_list$Des2)

#Now lets just grab the columns I need to make the array that PRAUC needs
PRAUC_Negative_True_L3 <- Noughty_list %>% select(P_CT_L3, Max_Read_count)
PRAUC_Negative_True_L2 <- Noughty_list %>% select(P_CT_L2, Max_Read_count)
PRAUC_Negative_True_L2 <- unique(PRAUC_Negative_True_L2)









### DESCARTES L2 ###

#Lets work on making the datatable I want- starting with the p_atlas_final_celltypes
#I want to translate the cell types from HPA into Descartes
#Lets read in the key CSV file
library(readr)
HPA_Des2_key <- read_csv("./Keys/HPA-Descartes_L2.csv")

#Now lets try to read and replace from the HPA Des3 Key
library(dplyr)
p_atlas_final_celltypes$Des2 <- ifelse(p_atlas_final_celltypes$combined %in% HPA_Des2_key$HPA, 
                                       HPA_Des2_key$`Descartes L2`[match(p_atlas_final_celltypes$combined, HPA_Des2_key$HPA)], 
                                       p_atlas_final_celltypes$combined[match(p_atlas_final_celltypes$combined, HPA_Des2_key$`Descartes L2`)])
p_atlas_final_celltypes <- p_atlas_final_celltypes[!duplicated(p_atlas_final_celltypes), ]

#Great - now we have a lovely table with the translated cell types for descartes L2.
#Now lets filter out all of the non 0 max values
Noughty_list_L2 <- p_atlas_final_celltypes %>% filter(Max_Read_count == '0')
#Lets just keep the proteins in my analysis
#Noughty_list <- p_atlas_final_celltypes %>% filter(Protein == c['FSHB','CABP2','FGF6',	'PTH',	'CEACAM18',	'IL3',	'OMP',	'IFNW1',	'MAGEA3',	'ALPI'])
#CT_protein_list <- read.csv(file=#path to your csv within quote marks "myfile.csv")
#For now lets just use: Positive_proteins
library(tidyverse)
Noughty_list_L2 <- filter(Noughty_list_L2,Protein %in% Positive_proteins$Var1)

#I need this to turn into an array with 2 columns I assume
#lets concatinate the protein and the Des3 cell type name
Noughty_list_L2$P_CT <- paste0(Noughty_list_L2$Protein,"_",Noughty_list_L2$Des2)

#Now lets just grab the columns I need to make the array that PRAUC needs
PRAUC_Negative_True <- Noughty_list_L2 %>% select(P_CT, Max_Read_count)














#NEGATIVE PREDICTION
#(Using the positive protein list)

# Load required libraries
library(dplyr)
library(readr)
library(tidyverse)

# Load protein and associated synapse keys
Synapse_key <- data.table::fread("./Synapse_keys.csv")
Synapse_key[, prot := sub("\\_.*", "", Name)]

# Load the list of proteins of interest
Positive_proteins <- read_csv("/Volumes/Expansion/R MAGMA/Pos_prot_2.csv")
#Positive_proteins <- read_csv("/Users/azamsaied/Documents/Research/Imperial/Nathan Skeene/Magma Celltyping/R MAGMA/Pos_prot_2.csv")

# Initialize the PM dataframe with the list of cell types
#PM <- read_csv("/Volumes/Expansion/R MAGMA/Protein_matrix.csv")
PM_L2 = matrix(, nrow = 0, ncol = 2)
PM_L3 = matrix(, nrow = 0, ncol = 2)
  
# Define the path and get list of files
#pth <- "/Volumes/Expansion/R MAGMA/CT_3000/Descartes L2/"
pth <- "/Volumes/Expansion/R MAGMA/CT_3000/Descartes L2/"
ss <- list.files(pth, pattern = "*.csv")
print(ss)
# Debugging print
print("Initial PM structure:")
print(head(PM_L2))

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
      temp <- temp %>% select(Celltype, P)
      temp$Celltype <- paste(prot_i,"_", temp$Celltype)
      temp <- temp %>% 
        group_by(Celltype) %>% 
        summarise(P = min(P))
      #temp$P <- (1-temp$P)
      temp$P <- (temp$P)
      
      # Rename the P column to the current protein name
      #temp <- temp %>% rename('Prediction' := 'P')
      temp <- temp %>% rename('P-Value' := 'P')
      
      # Debugging print
      print(paste("Merging data for protein:", prot_i))
      print(head(temp))
      
      # Merge PM and temp by Celltype
      #PM <- PM %>% left_join(temp, by = "Celltype")
      PM_L2 <- rbind(PM_L2, temp)
      
      # Debugging print to check the structure after merge
      print("PM structure after merge:")
      print(head(PM_L2))
      
    }, error = function(e) {
      print(paste0("Error in ", prot_i))
    })
  }
}

# Final structure check
print("Final PM structure:")
print(head(PM_L2))

# Save the final PM dataframe to a CSV file if needed
write_csv(PM_L2, "/Volumes/Expansion/R MAGMA/New Files/PRAUC_Negative_Prediction_L2.csv")
PRAUC_Negative_Prediction_L2 <- PM_L2
PRAUC_Negative_Prediction_L3 <- PM_L3
#And this bit gets rid of the spaces
PRAUC_Negative_Prediction_L2<- PRAUC_Negative_Prediction_L2 %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))

#And this bit gets rid of the spaces
PRAUC_Negative_Prediction_L3<- PRAUC_Negative_Prediction_L3 %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))

PRAUC_Negative_True_L2<- PRAUC_Negative_True_L2 %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))

PRAUC_Negative_True_L3<- PRAUC_Negative_True_L3 %>% 
  mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))

#Lets give them all the same column names
names(PRAUC_Negative_True_L2)[names(PRAUC_Negative_True_L2) == 'P_CT_L2'] <- 'Celltype'
names(PRAUC_Negative_True_L3)[names(PRAUC_Negative_True_L3) == 'P_CT_L3'] <- 'Celltype'

#Now we have to filter out the results which are not in our Negative True list

#PRAUC_Negative_Prediction_L3 <- PRAUC_Negative_Prediction_L3[row.names(PRAUC_Negative_Prediction_L3) %in% row.names(PRAUC_Negative_True_L3),]
#x_PRAUC_Negative_Prediction_L2 <- PRAUC_Negative_Prediction_L2[row.names(PRAUC_Negative_Prediction_L2) %in% row.names(PRAUC_Negative_True_L2),]
#dummy <- PRAUC_Negative_Prediction_L2$Celltype %in% PRAUC_Negative_True_L2$P_CT_L2
#dummy <- PRAUC_Negative_Prediction_L2[!PRAUC_Negative_Prediction_L2$Celltype %in% PRAUC_Negative_True_L2$Celltype,]


#XXXXXXXXXX-----CHECK THE FOLLOWING LINES OF CODE!!!!!----------XXXXXXXXXXX

x_PRAUC_Negative_Prediction_L2 <- PRAUC_Negative_Prediction_L2 %>% semi_join(PRAUC_Negative_True_L2, by = "Celltype")
x_PRAUC_Negative_Prediction_L3 <- PRAUC_Negative_Prediction_L3 %>% semi_join(PRAUC_Negative_True_L3, by = "Celltype")

write_csv(x_PRAUC_Negative_Prediction_L2, "/Volumes/Expansion/R MAGMA/PRAUC_Negative_Prediction_L2_filtered_0_proteins.csv")
write_csv(x_PRAUC_Negative_Prediction_L3, "/Volumes/Expansion/R MAGMA/PRAUC_Negative_Prediction_L3_filtered_0_proteins.csv")

head(dummy)


#Lets now get the positive predictions in there

PRAUC_Positive_Prediction_L2 <- read_csv("/Volumes/Expansion/R MAGMA/PRAUC_Positive_Prediction_L2.csv")
PRAUC_Positive_Prediction_L2 <- na.omit(PRAUC_Positive_Prediction_L2)
names(PRAUC_Positive_Prediction_L2)[names(PRAUC_Positive_Prediction_L2) == 'Cell_type-GWAS-L3'] <- 'Celltype'
names(PRAUC_Positive_Prediction_L2)[names(PRAUC_Positive_Prediction_L2) == 'P-value_L2'] <- 'P-Value'

PRAUC_Positive_Prediction_L3 <- read_csv("/Volumes/Expansion/R MAGMA/PRAUC_Positive_Prediction_L3.csv")
PRAUC_Positive_Prediction_L3 <- na.omit(PRAUC_Positive_Prediction_L3)
names(PRAUC_Positive_Prediction_L3)[names(PRAUC_Positive_Prediction_L3) == 'Cell_type-GWAS-L3'] <- 'Celltype'
names(PRAUC_Positive_Prediction_L3)[names(PRAUC_Positive_Prediction_L3) == 'P-value_L3'] <- 'P-Value'

library(data.table)

#To get PRROC to work  -I need to take -log of p-values
x_PRAUC_Negative_Prediction_L2$`P-Value` <- -log(x_PRAUC_Negative_Prediction_L2$`P-Value`)
x_PRAUC_Negative_Prediction_L3$`P-Value` <- -log(x_PRAUC_Negative_Prediction_L3$`P-Value`)
x_PRAUC_Positive_Prediction_L2 <- PRAUC_Positive_Prediction_L2
x_PRAUC_Positive_Prediction_L2$`P-Value` <- -log(x_PRAUC_Positive_Prediction_L2$`P-Value`)
x_PRAUC_Positive_Prediction_L3 <- PRAUC_Positive_Prediction_L3
x_PRAUC_Positive_Prediction_L3$`P-Value` <- -log(x_PRAUC_Positive_Prediction_L3$`P-Value`)

#OK - Lets finally do some Precision Recall analysis
#install.packages("PRROC")
library(PRROC)
packageVersion('PRROC')

roc_L2<-roc.curve(scores.class0 = x_PRAUC_Positive_Prediction_L2$`P-Value`, scores.class1 = x_PRAUC_Negative_Prediction_L2$`P-Value`, curve = TRUE)
pr_L2<-pr.curve(scores.class0 = x_PRAUC_Positive_Prediction_L2$`P-Value`, scores.class1 = x_PRAUC_Negative_Prediction_L2$`P-Value`, curve = TRUE)
roc_L2
plot(roc_L2)
pr_L2
plot(pr_L2)


roc_L3<-roc.curve(scores.class0 = x_PRAUC_Positive_Prediction_L3$`P-Value`, scores.class1 = x_PRAUC_Negative_Prediction_L3$`P-Value`, curve = TRUE)
pr_L3<-pr.curve(scores.class0 = x_PRAUC_Positive_Prediction_L3$`P-Value`, scores.class1 = x_PRAUC_Negative_Prediction_L3$`P-Value`, curve = TRUE)
roc_L3
plot(roc_L3)
pr_L3
plot(pr_L3)



#Because I've made all of these with p-values - I need to use a wrapper to get PRROC 
#to understand it'

#The wrapper works like this:
#pval_aucpr(pvals, causal_indexes, curve = FALSE)

#Step 1 - I need a vector of all of the p-values. I will put the true positives up top.
Pv_L2 <- rbind(PRAUC_Positive_Prediction_L2, x_PRAUC_Negative_Prediction_L2)
Pv_L3 <- rbind(PRAUC_Positive_Prediction_L3, x_PRAUC_Negative_Prediction_L3)

#Step 2 - Create the causal_index
causal_indexes <- 1:23

pval_aucpr(pvals, causal_indexes, curve = FALSE)


######    #######    #######    #####   ####   #####
#pr.curve( scores.class0, scores.class1=scores.class0, weights.class0=NULL,
#          weights.class1 = {if(is.null(weights.class0)){NULL}else{1-weights.class0}},
#          sorted = FALSE, curve = FALSE,
#          minStepSize=min(1,ifelse(is.null(weights.class0),1,sum(weights.class0)/100)),
#          max.compute=F, min.compute=F, rand.compute=F,dg.compute=T)

#scores.class0 the classification scores of i) all data points or ii) only the data points belonging
#to the positive class.
#scores.class1 the scores of the negative class if provided separately (see scores.class0)















































###Here is just a snippet of code to give me a manual list of proteins I need to 
#make sure we have Cell-typed with Descartes L3, because they are already on the L2 list

L2_pth <- "/Volumes/Expansion/R MAGMA/CT_3000/Descartes L2/"
L2_ss <- list.files(L2_pth, pattern = "*.csv")
L3_pth <- "/Volumes/Expansion/R MAGMA/CT_3000/Descartes L3/"
L3_ss <- list.files(L3_pth, pattern = "*.csv")

print(L2_ss)
print(L3_ss)

#Make an empty df
L2_prot_list<- c() 

for(ss_i in L2_ss) {
  print(paste("Processing file:", ss_i))
  ss_i_ne <- strsplit(strsplit(ss_i, ".csv")[[1]], "_")[[1]]
  prot_i <- ss_i_ne[[1]]
  syn_i <- ss_i_ne[[2]]
  #L2_prot_list <- bind_rows(L2_prot_list, prot_i)
  L2_prot_list <- rbind(L2_prot_list, prot_i)
}

#Make an empty df
L3_prot_list<- c() 

for(ss_i in L3_ss) {
  print(paste("Processing file:", ss_i))
  ss_i_ne <- strsplit(strsplit(ss_i, ".csv")[[1]], "_")[[1]]
  prot_i <- ss_i_ne[[1]]
  syn_i <- ss_i_ne[[2]]
  #L2_prot_list <- bind_rows(L2_prot_list, prot_i)
  L3_prot_list <- rbind(L3_prot_list, prot_i)
}

#write_csv(L3_prot_list, "/Volumes/Expansion/R MAGMA/L3_analysed_proteins.csv")
#write_csv(L2_prot_list, "/Volumes/Expansion/R MAGMA/L2_analysed_proteins.csv")


Prot_needed_to_celltype_with__L3 <- L2_prot_list[!(L2_prot_list %in% L3_prot_list),]
Prot_needed_to_celltype_with__L3 <- data.frame(Prot_needed_to_celltype_with__L3)
print(Prot_needed_to_celltype_with__L3)
write_csv(Prot_needed_to_celltype_with__L3, "/Volumes/Expansion/R MAGMA/Prot_needed_to_celltype_with_L3.csv")










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
