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
               "readr")

# set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs/1_Internal_validation")
if(!dir.exists(outpath)) dir.create(outpath)

# Step 1 Metabolomics Analytic Database Assembly ========================================
# Step 1.1 Import and Clean the Metabolomic Feature Table--------------------
# This workflow, we will focus on brain bank metabolome data that profiled recently in 2025yr

# Import the metabolome feature table: C18----------
mft <- read.delim(str_c(inpath, 
  "/Metabolomic_data/AIRCO_Liang_M625/c18neg",
  "/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt"),
  header = TRUE)

# View(mft %>% as_tibble())
# In this dataset, each row is a metabolomic feature-- hence, there are 10421 features in this dataset
dim(mft) #checks how many variables
# Starting from column #10, each column represents the signal intensity (relative concentrations)for 
mft %>% 
  dplyr::select(mz:Max.Intensity) %>% 
  head(5)
summary(mft[,1:10])
# filter out features with <15% detection rate and >30 median CV
mft %<>% 
  dplyr::filter(NumPres.All.Samples >= 0.15*ncol(mft-10) & median_CV <= 30)

# Import the metabolome feature table: HILIC---------
mft2 <- read.delim(str_c(inpath, 
  "/Metabolomic_data/AIRCO_Liang_M625/hilicpos",
  "/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt"),
  header = TRUE)
# View(mft2 %>% as_tibble())

# In this dataset, each row is a metabolomic feature-- hence, there are 13500 features in this dataset
dim(mft2) #checks how many variables
# Starting from column #10, each column represents the signal intensity (relative concentrations)for 
mft2 %>% 
  dplyr::select(mz:Max.Intensity) %>% 
  head(5)
summary(mft2[,1:10])
# filter out features with <15% detection rate and >30 median CV
mft2 %<>% 
  dplyr::filter(NumPres.All.Samples >= 0.15*ncol(mft2-10) & median_CV <= 30) # 2 features were filtered out


# Step 1.2 Extract Sample ID from mapping file-----------------------
# C18-------------------------------
# Import the sample sequence ID file
read.table(str_c(inpath, 
  "/Metabolomic_data/AIRCO_Liang_M625/c18neg",
  "/AIRCO_mapping_c18neg.txt"),
  header = TRUE, 
  fill = TRUE) %>% 
  set_names(c("File_Name", "Sample_ID", "Batch")) -> infofile
View(infofile)

# Extract file name as id.index
colnames(mft)[1:10]
colnames(mft)[11:dim(mft)[2]]
## create a vector of all of the column names 
colnames(mft)[11:dim(mft)[2]] %>% 
  as.data.frame() %>% set_names("File_Name") -> id.index
View(id.index)
# Extract File_Name and Sample_ID from id.index
infofile %>% 
  # as_tibble() %>% 
  dplyr::filter(File_Name %in% id.index$File_Name)  %>% 
  dplyr::mutate(Sample_ID = str_replace_all(Sample_ID, "_2", "")) %>% 
  dplyr::mutate(Sample_ID = case_when(str_detect(Sample_ID, "^nist|^q4") ~ str_c(Sample_ID,"_", Batch),
                TRUE ~ Sample_ID)) %>%
  dplyr::select(-Batch) -> var.list
# check duplicate
table(var.list$Sample_ID) %>% 
  as.data.frame() %>% 
  filter(Freq > 1) # none has two observations
                                
# make sure the file_names were correctly ordered with those in mft
var.list %<>% 
  dplyr::mutate(tmp = factor(File_Name, levels = id.index$File_Name)) %>% 
  arrange(tmp) %>% 
  dplyr::select(-tmp)
View(var.list)

# since IDs are totally numeric, we need to add a prefix to make sure they can be used as column names
var.list %<>% 
  dplyr::mutate(Sample_ID = case_when(!str_detect(Sample_ID, "^nist|^q4") ~ str_c("ID_", Sample_ID),
                                      TRUE ~ Sample_ID))

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

# HILIC-------------------------------
# Import the sample sequence ID file
read.table(str_c(inpath, 
  "/Metabolomic_data/AIRCO_Liang_M625/hilicpos",
  "/AIRCO_mapping_hilicpos.txt"),
  header = TRUE, 
  fill = TRUE) %>% 
  set_names(c("File_Name", "Sample_ID", "Batch")) -> infofile2
View(infofile2)

# Extract file name as id.index
colnames(mft2)[1:10]
colnames(mft2)[11:dim(mft2)[2]]
## create a vector of all of the column names 
colnames(mft2)[11:dim(mft2)[2]] %>% 
  as.data.frame() %>% set_names("File_Name") -> id.index2
View(id.index2)
# Extract File_Name and Sample_ID from id.index
infofile2 %>% 
  # as_tibble() %>% 
  dplyr::filter(File_Name %in% id.index2$File_Name)  %>% 
  dplyr::mutate(Sample_ID = str_replace_all(Sample_ID, "_1", "")) %>% 
  dplyr::mutate(Sample_ID = case_when(str_detect(Sample_ID, "^nist|^q4") ~ str_c(Sample_ID,"_", Batch),
                TRUE ~ Sample_ID)) %>%
  dplyr::select(-Batch) -> var.list2
# check duplicate
table(var.list2$Sample_ID) %>% 
  as.data.frame() %>% 
  filter(Freq > 1) # none has two observations
                                
# make sure the file_names were correctly ordered with those in mft
var.list2 %<>% 
  dplyr::mutate(tmp = factor(File_Name, levels = id.index2$File_Name)) %>% 
  arrange(tmp) %>% 
  dplyr::select(-tmp)
View(var.list2)

# since IDs are totally numeric, we need to add a prefix to make sure they can be used as column names
var.list2 %<>% 
  dplyr::mutate(Sample_ID = case_when(!str_detect(Sample_ID, "^nist|^q4") ~ str_c("ID_", Sample_ID),
                                      TRUE ~ Sample_ID))

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


# Step 1.3 Log transformation --------------------------------------
# C18 ----------------
dim(mdta)
###check metabolomics feature intensity distribution -- check the fifth metabolic feature
hist(as.numeric(mdta[30,11:459]))
hist(as.numeric(log2(mdta[30,11:459])))
###Apply universal log2 transformation to the feature table
mdta.log <- log2(mdta[,11:459])
mdta.log <- cbind(mdta[,1:10],mdta.log)
mdta.log %>% as_tibble()
# Convert all NaN values to NA
mdta.log %<>% 
  dplyr::mutate(across(everything(), ~ if_else(is.na(.x), NA, .x)))
View(mdta.log)

# HILIC ----------------
dim(mdta2)
###check metabolomics feature intensity distribution -- check the fifth metabolic feature
hist(as.numeric(mdta2[30,11:459]))
hist(as.numeric(log2(mdta2[30,11:459])))
###Apply universal log2 transformation to the feature table
mdta2.log <- log2(mdta2[,11:459])
mdta2.log <- cbind(mdta2[,1:10],mdta2.log)
mdta2.log %>% as_tibble()
# Convert all NaN values to NA
mdta2.log %<>% 
  dplyr::mutate(across(everything(), ~ if_else(is.na(.x), NA, .x)))
View(mdta2.log)

#  Upon this step, you have: ##################
## 1) filtered the feature table by removing low-presence and signals with bad data quality
## 2) replaced the colnames in the feature table with sample ID
## 3) replaced all missing values (0) with NA
## 4) log-transformed the data 
#  The feature table is clean and now ready for MWAS analysis
###############################################
# Step 1.4 Save and export datasets ----------------

# export the dataset
mdta.log %>% 
  write.csv(str_c(outpath, '/Met2025_mean.triplicate_log2.transformed_C18.csv'))
mdta2.log %>% 
  write.csv(str_c(outpath, '/Met2025_mean.triplicate_log2.transformed_HILIC.csv'))

# convert matrix from signature-per-row data into individual-per-row data ------
# C18---
mdta.log %>% 
  dplyr::select(11:ncol(mdta.log))  %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::mutate(Sample_ID = colnames(mdta.log[, 11:ncol(mdta.log)]), 
                GUID = case_when(str_detect(Sample_ID, "^ID_") ~ str_replace_all(Sample_ID, "ID_", ""), 
                                TRUE ~ Sample_ID)) %>% 
  dplyr::select(Sample_ID, GUID, everything()) -> m1
m1 %>% 
  write.csv(str_c(outpath, '/Met2025_C18_for.MWAS.csv'))

# HILIC---
mdta2.log %>% 
  dplyr::select(11:ncol(mdta2.log))  %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::mutate(Sample_ID = colnames(mdta2.log[, 11:ncol(mdta2.log)]), 
                GUID = case_when(str_detect(Sample_ID, "^ID_") ~ str_replace_all(Sample_ID, "ID_", ""), 
                                TRUE ~ Sample_ID)) %>% 
  dplyr::select(Sample_ID, GUID, everything()) -> m2
m2 %>% 
  write.csv(str_c(outpath, '/Met2025_HILIC_for.MWAS.csv'))


# Step 2: Map metabolites to be validated ######################################################
library(data.table)
# match metabolites to be validated ----------------
mdta.log = fread(str_c(outpath, "/Met2025_mean.triplicate_log2.transformed_C18.csv"))
mdta2.log = fread(str_c(outpath, "/Met2025_mean.triplicate_log2.transformed_HILIC.csv"))

# import metabolite to be validated ---------
meta.index = read_xlsx(str_c(work.dir, '/0_Manuscript/Send_to_coauthors/2_Supplementary Data Tables_Apr01.2026.xlsx'), sheet = "Table S5_identified metabolites")
meta.index %>% 
  dplyr::select(all_of(c("m/z", "rt (s)", "Mode", "Metabolite name"))) %>% 
  set_names(c("mz", "time", "Mode", "Metabolite_name")) -> meta.index

# import orginal variable names (as pair index)
tmp =  readxl::read_xlsx(str_c(work.dir, "/0_Manuscript/tables/PN_Identified metabolites (Level-1 evidence).xlsx")) %>% 
  dplyr::rename(mz = `m/z`, time = `rt (second)`) %>% 
  dplyr::select(mz, time, Mode, Components)

meta.index %<>% 
  left_join(tmp, by = c("mz", "time", "Mode"))

View(meta.index)
dim(mdta.log)
dim(mdta2.log)
meta_match = function(i){
  if(meta.index$Mode[i]=="Negative"){
    # c18
    mdta.log %>% as_tibble() %>% 
      dplyr::filter(mz >= (meta.index$mz[i] - (meta.index$mz[i]/1000000*10)) & mz <= (meta.index$mz[i] + (meta.index$mz[i]/1000000*10))) %>% # 10 ppm
      #dplyr::filter(time >= (meta.index$time[i]-60) & time <= (meta.index$time[i]+60)) %>% 
      dplyr::mutate(mode = "Negative", 
                    metabolite = meta.index$Metabolite_name[i], 
                    Components = meta.index$Components[i]) %>% 
      dplyr::select(mode, metabolite, Components, everything()) -> tmp
  }else if(meta.index$Mode[i]=="Positive"){
    # hilic
    mdta2.log %>% as_tibble() %>% 
      dplyr::filter(mz >= (meta.index$mz[i] - (meta.index$mz[i]/1000000*10)) & mz <= (meta.index$mz[i] + (meta.index$mz[i]/1000000*10))) %>% # 10 ppm
      #dplyr::filter(time >= (meta.index$time[i]-60) & time <= (meta.index$time[i]+60)) %>% 
      dplyr::mutate(mode = "Positive", 
                    metabolite = meta.index$Metabolite_name[i], 
                    Components = meta.index$Components[i]) %>% 
      dplyr::select(mode, metabolite, Components, everything()) -> tmp
  }
  return(tmp)
}

nrow(meta.index) # 19 metabolic features in total to be matched

# 10ppm
# no retention time threshold set
match = purrr::map_dfr(1:nrow(meta.index), meta_match) # 18 matches, where multiple matching exist

match  %>% colnames()
# delete multiple matching according to meta.index
# we will delete row=1 and row=13 since they are duplicates
match[-c(1,13), ] %>% 
  write.csv(str_c(outpath, '/1_Met2025_metabolites_to_be_modelled.csv'))


