rm(list=(ls()))
# 1- Prepare steps =====================================
# load package ---------
pacman::p_load("readxl", 
               "purrr", 
               "lubridate", 
               "dplyr", 
               "magrittr",
               "stringr", 
               "rlist", 
               "tidyverse", 
               "lmerTest", 
               "lme4", 
               "readr", 
               "data.table")
# set work directory----------
getwd() # returns "C:/Users/MREN30/OneDrive - Emory/1_Ongoing_project/2_TRAP_ADRD"
# # so, you need to set work directory as something like 
# setwd("C:/.../2_TRAP_ADRD" )
setwd("C:/Users/MREN30/OneDrive - Emory/1_Ongoing_project/2_TRAP_ADRD")
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# Step 1 Metabolomics Analytic Database Assembly ========================================
# Step 1.1 Import and Clean the Metabolomic Feature Table--------------------
# Import the metabolome feature table: C18----------
mft <- fread(str_c(inpath, 
  "/Metabolomic_data",
  "/c18neg.csv"),
  header = TRUE)
# have direct view of data
# View(mft %>% as_tibble())
# In this dataset, each row is a metabolomic feature-- hence, there are 15927 features in this dataset
dim(mft) #checks how many variables(rows) and samples(~= columns) of this dataset
# Starting from column #11, each column represents the signal intensity (relative concentrations) of one individual participant
#have a look at the information on metabolomic signatures of negative mode (c18)
mft %>% 
  dplyr::select(mz:Max.Intensity) %>% 
  head(5)
summary(mft[,1:10])
# To guarantee we use metabolomic sigatures with high quality 
# (rather than signatures "not real biological signatures, but the noises that came from LCMS"), we need to
# filter out features with <15% detection rate and >30 median coefficient variation (CV)
mft %<>% 
  dplyr::filter(NumPres.All.Samples >= 0.15*ncol(mft-10) & median_CV <= 30) 
nrow(mft) # 2 features were removed out

# Import the metabolome feature table: HILIC---------
mft2 <- fread(str_c(inpath, 
  "/Metabolomic_data",
  "/hilicpos.csv"),
  header = TRUE)
# View(mft2 %>% as_tibble())
# In this dataset, each row is a metabolomic feature-- hence, there are 20,051 features in this dataset
dim(mft2) #check rows and columns (20051 versus 197)
# have a look at signatures information of positive mode (hilic)
mft2 %>% 
  dplyr::select(mz:Max.Intensity) %>% 
  head(5)
summary(mft2[,1:10])
# filter out features with <15% detection rate and >30 median CV
mft2 %<>% 
  dplyr::filter(NumPres.All.Samples >= 0.15*ncol(mft2-10) & median_CV <= 30)
nrow(mft2) # 10 features were filtered out

# Step 1.2 Extract Sample ID from mapping file-----------------------
# 1.2.1 negative mode: C18-------------------------------
# Import the sample sequence ID file
read.table(str_c(inpath, 
  "/Metabolomic_data",
  "/ADRC_mapping_c18neg.txt"),
  header = TRUE, 
  fill = TRUE) %>% 
  set_names(c("File_Name", "Sample_ID", "Batch")) -> infofile
# View(infofile)
# Extract file name as id.index
colnames(mft)[1:10] #just a quick look
colnames(mft)[11:dim(mft)[2]] #just a quick look
## create a vector of all of the column names called in.index
colnames(mft)[11:dim(mft)[2]] %>% 
  as.data.frame() %>% set_names("File_Name") -> id.index
# View(id.index)
# Extract File_Name and Sample_ID from id.index
infofile %>% 
  dplyr::filter(File_Name %in% id.index$File_Name)  %>% 
  dplyr::mutate(Sample_ID = str_extract(Sample_ID, "(\\w(\\d+|\\w+)-\\d+)|nist\\d|(q3\\w+\\d+_(\\d+\\w))"), 
                Sample_ID = str_replace(Sample_ID, "-", "_")) %>% 
  dplyr::select(-Batch) -> var.list

# check duplicate
table(var.list$Sample_ID) %>% 
  as.data.frame() %>% 
  filter(Freq > 1)
# arrange var.list dataset, to make sure the file_names were correctly ordered with those in mft
var.list %<>% 
  dplyr::mutate(tmp = factor(File_Name, levels = id.index$File_Name)) %>% 
  arrange(tmp) %>% 
  dplyr::select(-tmp)
# View(var.list)
# rename colnames of mft
mft %>% 
  dplyr::select(any_of(id.index$File_Name)) %>% 
  set_names(var.list$Sample_ID)  %>% as_tibble() %>% 
  dplyr::mutate(across(everything(), ~dplyr::na_if(.x, 0))) -> tmp1
# combine metabolomic signature information
mft %>%
  as_tibble() %>% 
  dplyr::select(-any_of(id.index$File_Name)) %>% 
  cbind(tmp1) -> mdta
# View(mdta)

# 1.2.2 positive mode: HILIC----------------------------------------------
# Import the sample sequence ID file
read.table(str_c(inpath, 
  "/Metabolomic_data",
  "/ADRC_mapping_hilicpos.txt"),
  header = TRUE, 
  fill = TRUE) %>% 
  set_names(c("File_Name", "Sample_ID", "Batch")) -> infofile2
# View(infofile2)
# Extract file name as id.index
colnames(mft2)[1:10]
colnames(mft2)[11:dim(mft2)[2]]
## create a vector of all of the column names 
colnames(mft2)[11:dim(mft2)[2]] %>% 
  as.data.frame() %>% set_names("File_Name") -> id.index2
# View(id.index2)
# Extract File_Name and Sample_ID from id.index
infofile2 %>% 
  dplyr::filter(File_Name %in% id.index2$File_Name)  %>% 
  dplyr::mutate(Sample_ID = str_extract(Sample_ID, "(\\w(\\d+|\\w+)-\\d+)|nist\\d|(q3\\w+\\d+_(\\d+\\w))"), 
                Sample_ID = str_replace(Sample_ID, "-", "_")) %>% 
  dplyr::select(-Batch) -> var.list2
# check duplicates
table(var.list2$Sample_ID) %>% 
  as.data.frame() %>% 
  filter(Freq > 1)# E0586 has two observations
# make sure the file_names were correctly ordered with those in mft
var.list2 %<>% 
  dplyr::mutate(tmp = factor(File_Name, levels = id.index2$File_Name)) %>% 
  arrange(tmp) %>% 
  dplyr::select(-tmp)
# View(var.list2)
# rename colnames of mft
mft2 %>% 
  dplyr::select(any_of(id.index2$File_Name)) %>% 
  set_names(var.list2$Sample_ID)  %>% as_tibble() %>% 
  dplyr::mutate(across(everything(), ~dplyr::na_if(.x, 0))) -> tmp2
# combine metabolomic signature information
mft2 %>%
  as_tibble() %>% 
  dplyr::select(-any_of(id.index2$File_Name)) %>% 
  cbind(tmp2) -> mdta2
# View(mft2)

# Step 1.3 Log transformation --------------------------------------
# 1.3.1 C18 ----------------
dim(mdta) # 15925 198
# check metabolomics feature intensity distribution 
# check the 20th metabolic feature as an example
hist(as.numeric(mdta[20,11:dim(mdta)[2]]))
hist(as.numeric(log2(mdta[20,11:dim(mdta)[2]])))
###Apply universal log2 transformation to the feature table
mdta.log <- log2(mdta[,11:dim(mdta)[2]])
mdta.log <- cbind(mdta[,1:10],mdta.log)
# Convert all NaN values to NA
mdta.log %<>% 
  dplyr::mutate(across(everything(), ~ if_else(is.na(.x), NA, .x)))
# View(mdta.log)

# 1.3.2 HILIC ----------------
dim(mdta2)
###check metabolomics feature intensity distribution -- check the fifth metabolic feature
hist(as.numeric(mdta2[20,11:dim(mdta2)[2]]))
hist(as.numeric(log2(mdta2[20,11:dim(mdta2)[2]])))
###Apply universal log2 transformation to the feature table
mdta2.log <- log2(mdta2[,11:dim(mdta2)[2]])
mdta2.log <- cbind(mdta2[,1:10],mdta2.log)
# Convert all NaN values to NA
mdta2.log %<>% 
  dplyr::mutate(across(everything(), ~ if_else(is.na(.x), NA, .x)))
# View(mdta2.log)

#  Upon this step, you have: ##################
## 1) filtered the feature table by removing low-presence and signals with bad data quality
## 2) replaced the colnames in the feature table with sample ID
## 3) replaced all missing values (0) with NA
## 4) log-transformed the data 
#  The feature table is clean and now ready for MWAS analysis
###############################################

# Step 1.5 Save and export datasets ----------------
# export the dataset
mdta.log %>% 
  fwrite(str_c(outpath, '/log2.transformed.metabolomics_C18.csv'))
mdta2.log %>% 
  fwrite(str_c(outpath, '/log2.transformed.metabolomics_HILIC.csv'))

# convert matrix from signature-level data into individual-level data ------
# C18--
mdta.log %>% 
  dplyr::select(11:ncol(mdta.log))  %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::mutate(Sample_ID = colnames(mdta.log[, 11:ncol(mdta.log)])) %>% 
  dplyr::select(Sample_ID, everything()) %>% 
  fwrite(str_c(outpath, '/Cleaned_data_for_MWAS_c18.csv'))
# HILIC--
mdta2.log %>% 
  dplyr::select(11:ncol(mdta2.log))  %>%
  t() %>% 
  as_tibble() %>% 
  dplyr::mutate(Sample_ID = colnames(mdta2.log[, 11:ncol(mdta2.log)])) %>% 
  dplyr::select(Sample_ID,  everything())  %>% 
  fwrite(str_c(outpath, '/Cleaned_data_for_MWAS_hilic.csv'))

# save work environment as rdata file. 
save.image(str_c(getwd(), "/3_Outputs",
                          "/!Data_for_data_processing_3.rdata"))

# # next time when you want to retrieve any results that calculated with this code file, just run code load() to call the environment
# load(str_c(getwd(), "/3_Outputs",
#                     "/!Data_for_data_processing_3.rdata"))

