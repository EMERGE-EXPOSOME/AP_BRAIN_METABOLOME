rm(list=(ls()))
# 1- Prepare steps =====================================
# load package ---------
pacman::p_load("readxl", 
               "purrr",
               "future",
               "lubridate", 
               "dplyr", 
               "magrittr",
               "stringr", 
               "rlist", 
               "tidyverse", 
               "data.table")

# set work directory-----------------------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# 2- Import demographic data-----------------------
fread(str_c(inpath, "/Demographics/Brain_bank_Marla_clean_432obs.csv")) -> df
View(df)
# 2.1 data clean for covaraites and outcome variables------------------------

### 2.1.1 Age [need discussion]----
# age at death
df$Age_at_Death_Bx
df$age_at_Death_new %>% summary()
df %<>% 
  dplyr::mutate(Age_MY = case_when(is.na(age_at_Death_new) & !is.na(Age_at_Death_Bx) ~ as.numeric(str_extract(Age_at_Death_Bx, "\\d+")), 
                                   TRUE ~ age_at_Death_new))
# check the age distribution, Age_MY: Age variable used by Mengyuan and Yasmith
# We can see that which is not reasonable for the age of death
df$Age_MY %>% summary() 

# age onset
df$age_at_onset_new %>% summary()
# We can see that there remains a problem, which is: there are some small age values with a minimum value of 10
# which is not reasonable for the age of onset
df %>% 
  dplyr::select(Duration_.years., duration_cal) %>% View()
df$Duration_.years.
df$duration_cal %>% summary()
# in line 189, there is a discrepancy between the duration_cal and Duration_.years.
# which we should use as final value?

### 2.1.2 BMI [need discussion] ----
# we don't have weight or BMI data for this population??
# Answer is yes, we don't have BMI information

### 2.1.3 Sex  ----
df$sex %>% table(useNA = "always") 
df$sex = as.factor(df$sex)

### 2.1.4 race  ----
df$race %>% table(useNA = "always")
df$race = as.factor(df$race)

### 2.1.5 education ----
df$edu_cat %>% table(useNA = "always")
df %<>% 
  dplyr::mutate(across(edu_cat, ~ factor(as.factor(.x), levels = c("High school or less", "Graduate degree", "College degree", NA))))

### 2.1.6 Family history of AD ----
df$family_hx_new %>% table(useNA = "always")
df$family_hx_new = as.factor(df$family_hx_new)

### 2.1.7 PMI(Postmortem interval) [need discussion]----
df$PMI_.hr. %>% table(useNA = "always")
df$PMI_hr_num %>% summary()
df %>% 
  dplyr::select(PMI_.hr., PMI_hr_num) %>% 
  dplyr::filter(is.na(PMI_hr_num) & !str_detect(PMI_.hr., "na|minutes")) %>% View()
# you can see that there are still some values like "prob<5", prob ~12, were not converted to proxy numerics
# instead, we can
df %<>% 
  dplyr::mutate(PMI_MY = case_when(is.na(PMI_hr_num) & !str_detect(PMI_.hr., "na|minutes") ~ as.numeric(str_extract(PMI_.hr., "\\d+")), 
                                   TRUE ~ PMI_hr_num))
df$PMI_MY %>% summary()
df %>% 
  dplyr::select(PMI_.hr., PMI_hr_num, PMI_MY) %>% 
  dplyr::filter(is.na(PMI_hr_num) & !str_detect(PMI_.hr., "na|minutes")) %>% View()
### 2.1.8 ADI(Area Deprivation Index) ----
df$ADI_NATRANK %>% summary()

### 2.1.9 APOE ----
df %>% 
  dplyr::select(ApoE, ApoE_bi, ApoE_num, ApoE_tri, Case_Number_APOE, CRIN_APOE, APOE_FINAL) %>% 
  dplyr::filter(Case_Number_APOE == "NA") %>% View()
df$APOE_FINAL %>% table(useNA = "always") # specific ApoE status

df %<>% 
  dplyr::mutate(across(ApoE_bi, ~ factor(as.factor(.x), labels = c("E4 absence", "E4 presence"))))

df$ApoE_bi %>% table(useNA = "always") # binary ApoE status including E4 absence or presence


df %<>% 
  dplyr::mutate(across(ApoE_tri, ~ factor(as.factor(.x), labels = c("E4 absence", "E4 heterozygous", "E4 homozygous"))))
df$ApoE_tri %>% table(useNA = "always") # tri-class ApoE status including 0=E4 absence, 1=E4 heterozygous, 2=E4 homozygous 

df %>% 
  dplyr::select(ApoE, ApoE_bi, ApoE_num, ApoE_tri, Case_Number_APOE, CRIN_APOE, APOE_FINAL) 


### 2.1.10 ADRD related outcomes [need discussion] ----
# ABC (Alzheimer's Biomarker Consortium) score
df %<>% 
  dplyr::mutate(across(ABC_new, ~ factor(as.factor(.x), levels = c("Not", "Low", "Intermediate", "High", NA))))
df$ABC_new %>% table(useNA = "always") # ABC categorical variable
df$ABC_num %>% table(useNA = "always") # ABC numerical variable

# Braak stage (Braak stage of AD pathology)[need discussion]
df$Braak_Stage %>% table(useNA = "always")
df$Braak_Stage_new %>% table(useNA = "always") # Braak stage categorical variable, is there any way to combine these stages in clinical meaningful perspective?
df$Braak_Stage_num %>% table(useNA = "always") # Braak stage numerical variable

# CERAD (Consortium to Establish a Registry for Alzheimer's Disease) score
df$CERAD %>% table(useNA = "always")
df %<>% 
  dplyr::mutate(across(CERAD_new, ~ factor(as.factor(.x), levels = c("No", "Sparse", "Moderate", "Frequent", NA))))
df$CERAD_new %>% table(useNA = "always")
df$CERAD_num %>% table(useNA = "always")

# others: Thal, ADReagan [need discussion]
# Thal (Thal phase of AD pathology)
df$Thal %>% table(useNA = "always") 
df$Thal_new %>% table(useNA = "always") # Thal stage categorical variable, is there any way to combine some of these stages in clinical meaningful perspective?

# ADReagan (Reagan stage of AD pathology)
df$ADReagan %>% table(useNA = "always")
df %<>% 
  dplyr::mutate(across(ADReagan_new, ~ factor(as.factor(.x), levels = c("Not", "Low", "Intermediate", "High", NA))))
df$ADReagan_new %>% table(useNA = "always")

### 2.1.11 Participant IDs variables ----
df$Case_Number
df$OBJECTID
df$GUID
df$CRIN_number

### 2.1.12 frozen tissue ---- 
df$Frozen_Tissue_new %>% table(useNA = "always") # X represents available frozen tissue. 

### 2.1.13 Neuropathology data ----
df$Clin_Dx %>% table(useNA = "always")
df$Primary_Neuropathologic_Diagnosis %>% table(useNA = "always")
df$Secondary_Neuropathologic_Diagnosis %>% table(useNA = "always")
df$Tertiary_Neuropathologic_Diagnosis %>% table(useNA = "always")
df$dement_class %>% table(useNA = "always")

# 2.2 Export cleaned demographc table------------------------
df %>% 
  dplyr::select(
    Case_Number, OBJECTID, GUID, CRIN_number, # subject id and grid id
    Age_MY,age_at_Death_new, Age_at_Death_Bx, # age at death
    age_at_onset_new, # age at onset
    sex, race, edu_cat, family_hx_new, # demographic variables
    PMI_MY, # PMI, suppose that PMI_MY is the final PMI variable
    ADI_NATRANK, # ADI
    APOE_FINAL, ApoE_bi, ApoE_tri, # ApoE status
    ABC_new, ABC_num, # ABC score
    Braak_Stage_new, Braak_Stage_num, # Braak stage
    CERAD_new, CERAD_num, # CERAD score
    Thal_new, # Thal phase
    ADReagan_new, # ADReagan stage
    Frozen_Tissue_new, # frozen tissue
    dement_class, Clin_Dx, Primary_Neuropathologic_Diagnosis, Secondary_Neuropathologic_Diagnosis, Tertiary_Neuropathologic_Diagnosis # Neuropathology data
  ) -> exp
exp %>% 
  fwrite(str_c(outpath, "/Cleaned_data_demographics.csv"))






 