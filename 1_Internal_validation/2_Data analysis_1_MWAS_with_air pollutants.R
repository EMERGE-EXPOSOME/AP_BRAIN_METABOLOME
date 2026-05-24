rm(list=(ls()))
# 1- Prepare steps -------------------------------------
# load package ---------
pacman::p_load("readxl", 
               "purrr", 
               "lubridate", 
               "dplyr", 
               "magrittr",
               "stringr", 
               "rlist", 
               "tidyverse", 
               "tableone", 
               "gghalves", 
               "ddpcr",
               "ggpubr", 
               "corrplot", 
               "glue", 
               "qgcomp", 
               "data.table", 
               "lme4", 
               "lmerTest", 
               "statmod",
               "furrr", 
               "future")
# set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs/1_Internal_validation")
# outpath <- str_c(getwd(), "")

# Create output folder for the upcoming MWAS results of traffic-related air pollution and metabolomic signatures
if(!file.exists(str_c(outpath,"/AP_MWAS"))){
  dir.create(str_c(outpath,"/AP_MWAS"))
}

# 2-Load data ----------
# load demographic data
fread(str_c(work.dir, "/3_Outputs/Cleaned_data_demographics.csv")) -> demo
demo %<>% 
  dplyr::mutate(ApoE_tri = if_else(ApoE_tri == "", NA, ApoE_tri))
demo$GUID

# load air pollution data
fread(str_c(work.dir, "/3_Outputs/Cleaned_data_Air_pollution.csv")) -> ap.data
ap.data$GUID

# load metabolomic data
fread(str_c(outpath, "/1_Met2025_metabolites_to_be_modelled.csv")) %>%
  dplyr::select(-V1) %>% 
  as_tibble() %>%
  dplyr::mutate(V1 = row_number()) -> m.data

dim(m.data)

m.data %>% 
  dplyr::select(14:ncol(m.data))  %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::mutate(Sample_ID = colnames(m.data[, 14:ncol(m.data)]), 
                GUID = case_when(str_detect(Sample_ID, "^ID_") ~ str_replace_all(Sample_ID, "ID_", ""), 
                                TRUE ~ Sample_ID)) %>% 
  dplyr::select(GUID, everything(), -Sample_ID) %>% 
  dplyr::filter(!str_detect(GUID, "^q4|^nist")) -> meta.data 

dim(meta.data)
# impute missings as mim/2
meta.data %<>%
  dplyr::mutate(across(all_of(matches("V\\d+")), ~ifelse(is.na(.x), min(.x/2, na.rm = TRUE), .x)))

meta.data$GUID = as.integer(meta.data$GUID)

colnames(demo)
colnames(meta.data)


# merge data
demo %>% 
  left_join(ap.data, by = "Case_Number") -> tmp
colnames(tmp)

tmp %>% 
  dplyr::rename(GUID = "GUID.x") %>%
  left_join(meta.data, by = "GUID") %>% 
  as_tibble() -> df.vali

# save data
save.image(str_c(outpath, "/!Data_for_validation.rdata"))
###You are now ready for the MWAS data analysis----
load(str_c(outpath, "/!Data_for_validation.rdata"))


df.vali %>% nrow()
colnames(df.vali)[1:150]


# 
df.vali %>% 
  dplyr::filter(!is.na(V1)) %>% # those have metabolomics data -- 347 in total
  dplyr::filter(!is.na(.[["SEDAC_PM25NO3_1y"]])) -> tmp # those have air pollution data -- 300 in total

# import those have been analyzed in discovery phase
df.dis = fread(str_c(work.dir, "/3_Outputs/Cleaned_TRAP_data_for_brain_donor_N162.csv"))
# GUID in discovery cohort
df.dis %>% 
  dplyr::filter(!is.na(SEDAC_PM25NO3_1y)) %>% 
  .$GUID -> guid.dis

# final data for validation
tmp %>% 
  dplyr::filter(!(GUID %in% guid.dis)) -> df.validation #143 obs for validation

colnames(df.validation)


###calculate demographics######
tmp %>% 
  dplyr::mutate(flag = case_when(GUID %in% guid.dis ~ "Discovery cohort", 
                                 TRUE ~ "Validation cohort")) -> tmp.2

# of note: in 2.3.1, here we did not exclude participants without TRAP data
tt = c("Age_MY", "sex", "race", "edu_cat", "PMI_MY",
       "ADI_NATRANK", "ApoE_tri", "ABC_new", 
       "Braak_Stage_new", "CERAD_new") 
cat = c("sex", "race", "edu_cat", 
       "ApoE_tri", "Braak_Stage_new") 

# create TableOne: AllPTBcsAllTerm
tab1 <-  print(CreateTableOne(vars = tt, # here we use All term birth as control group, to show whole participants in our study
                              data = tmp.2, factorVars = cat,
                              includeNA = FALSE, 
                              strata = "flag"),
               showAllLevels = TRUE)
rownames(tab1)  %>% as_tibble()  %>%
  cbind(tab1 %>% as_tibble()) %>%
  as_tibble() %>%
  set_names(c("Variables", "Levels","Discovery cohort", "Validation cohort", "p-value", "test")) %>% 
writexl::write_xlsx(stringr::str_c(outpath,
                                  "/TableOne_1_Descriptive_results.xlsx"))

tmp.2$edu_cat %>% table(useNA = "ifany")

##############################




# 3- load pairs that need to be verified in validation analyses ----------------------
m.data %<>% 
  dplyr::mutate(New.name = str_c("V", as.character(row_number()))) %>% 
  dplyr::select(mode, metabolite, Components, New.name, everything())# m.data is used to link 
m.data %<>% 
  dplyr::rename(Mode = "mode")

# import metabolite to be validated ---------
fread(str_c(work.dir, "/3_Outputs/DataAnalysis_6_Sensitivity_analysis/SEDAC_individual_pollutant_sensitivitiy_analysis_Adding_ABC_new.csv")) %>% 
  dplyr::filter(Group == "Main analysis") %>% 
  dplyr::select(Exposure, Components, Mode) %>% # all were significant associations, that is why they were in sensitivity analysis group
  dplyr::mutate(Exposure = str_c("SEDAC_", Exposure)) %>% 
  left_join(m.data %>% dplyr::select(Mode, metabolite, Components, New.name), by = c("Components", "Mode")) -> pairs

pairs %<>% 
  dplyr::filter(!is.na(New.name))


# 4- Metabolomics-wide Association Study (MWAS)------------------------------
i=1

Ex.vali = function(i){
  
  df = data %>% 
    dplyr::filter(!is.na(.[[pairs$Exposure[i]]]))

  print(glue("For TRAP [ {pairs$Exposure[i]} ], a total of [ {nrow(df)} ] obs were obtained for downstream analysis!"))
  

  IQR = quantile(as.matrix(df[,pairs$Exposure[i]]),0.75) - quantile(as.matrix(df[,pairs$Exposure[i]]),0.25)
  # adjusted aCoef and 95% CI
  eqt = paste(pairs$New.name[i], " ~ ", pairs$Exposure[i], " + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri")

  fit <- lm(as.formula(eqt), data = df)
  tmp <- summary(fit)$coefficients
  if(is.null(tmp) | !(pairs$Exposure[i] %in% row.names(tmp))){
    # in some cases, model may not fit (due to low variance of specific metabolites)
    # so we need to assign tmp1 as follow:
    tmp2 <- cbind(aCoef = NA,
                  LCI_a = NA,
                  UCI_a = NA, 
                  t_a = NA,
                  Pvalue_a = NA)
    print(str_c("Warning!! metabolite [", pairs$metabolite[i], "] could not fit adjusted models for ", pairs$Exposure[i], " due to low variates !"))
  }else{      
    tmp2 <- cbind(aCoef = coef(fit)[[2]]*IQR,
                  LCI_a = confint.default(fit)[2,1]*IQR,
                  UCI_a = confint.default(fit)[2,2]*IQR, 
                  t_a = tmp[2,3],
                  Pvalue_a = tmp[2,4])
  }

  #3. combine cCoef and aCoef
  tmp2 %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(Exposure = pairs$Exposure[i],
                  Metabolite = pairs$metabolite[i], 
                  N_vali = nrow(df)) %>% 
    dplyr::select(Exposure, Metabolite, N_vali, everything()) -> res

  return(res)
}


# execute
data = df.validation
tictoc::tic()
res = purrr::map_dfr(1:nrow(pairs), Ex.vali)
tictoc::toc()

res %>% 
  fwrite(str_c(outpath, "/AP_MWAS/1_complete_new_data(N143).csv"))




# execute total sample size
data = tmp
tictoc::tic()
res2 = purrr::map_dfr(1:nrow(pairs), Ex.vali)
tictoc::toc()
res2 %>% View()

res2 %>% 
  fwrite(str_c(outpath, "/AP_MWAS/1_all_data(N300).csv"))
