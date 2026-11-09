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
               "statmod"
               )

# 1- set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# load data
load(str_c(outpath, "/!Data_for_Data_Analysis_1.rdata"))
ls()

# read metabolite list to be used
metabolites <- fread(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_fdr0.2_final.csv")) %>% 
              dplyr::select(mz, time, Components, Mode, metabolite_1)

# read pairs need to be modeled for stratification analyses
pairs <- fread(str_c(outpath, "/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
                                 "/Annotated_metabolites_XtoM_SEDAC_thresholds_mz10ppm_fdr0.2.csv"))

# 2- features used for stratification --------------------
df.c18$Age_MY
df.c18$sex %>% table()
df.c18$ApoE_tri %>% table(useNA = "always")
df.c18$ABC_new %>% table(useNA = "always")
df.c18$Braak_Stage_new %>% table(useNA = "always")
df.c18$CERAD_new %>% table(useNA = "always")

df.c18 %<>% 
  dplyr::mutate(
                Age_cat65 = case_when(Age_MY < 65 ~ "<65 yr", 
                                    Age_MY >= 65 ~ ">=65 yr", 
                                    TRUE ~ NA),     
                Age_cat75 = case_when(Age_MY < 75 ~ "<75 yr", 
                                    Age_MY >= 75 ~ ">=75 yr", 
                                    TRUE ~ NA),     
                Age_cat = case_when(Age_MY < 70 ~ "<70 yr", 
                                    Age_MY >= 70 ~ ">=70 yr", 
                                    TRUE ~ NA), 
                ApoE_cat = case_when(ApoE_tri  %in% c("E4 homozygous", "E4 heterozygous") ~ "E4 existance", 
                                     ApoE_tri == "E4 absence" ~ "E4 absence", 
                                     TRUE ~ NA), 
                ABC_cat = case_when(ABC_new %in% c("High", "Intermediate") ~ "High", 
                                    ABC_new %in% c("Low", "Not") ~ "Low",
                                    ABC_new == "" ~ NA, 
                                    TRUE ~ NA), 
                Braak_cat = case_when(Braak_Stage_new %in% c("Stage 0", "Stage 1", "Stage 2", "Stage 3") ~ "Lower stage", 
                                      Braak_Stage_new %in% c("Stage 4", "Stage 5", "Stage 6") ~ "Higher stage", 
                                      Braak_Stage_new == "" ~ NA, 
                                      TRUE ~ NA), 
                CERAD_cat = case_when(CERAD_new %in% c("Frequent", "Moderate") ~ "Frequent", 
                                      CERAD_new %in% c("No", "Sparse") ~ "Sparse", 
                                      TRUE ~ NA))
df.hilic %<>% 
  dplyr::mutate(
                Age_cat65 = case_when(Age_MY < 65 ~ "<65 yr", 
                                    Age_MY >= 65 ~ ">=65 yr", 
                                    TRUE ~ NA),     
                Age_cat75 = case_when(Age_MY < 75 ~ "<75 yr", 
                                    Age_MY >= 75 ~ ">=75 yr", 
                                    TRUE ~ NA), 
                Age_cat = case_when(Age_MY < 70 ~ "<70 yr", 
                                    Age_MY >= 70 ~ ">=70 yr", 
                                    TRUE ~ NA), 
                ApoE_cat = case_when(ApoE_tri  %in% c("E4 homozygous", "E4 heterozygous") ~ "E4 existance", 
                                     ApoE_tri == "E4 absence" ~ "E4 absence", 
                                     TRUE ~ NA), 
                ABC_cat = case_when(ABC_new %in% c("High", "Intermediate") ~ "High", 
                                    ABC_new %in% c("Low", "Not") ~ "Low",
                                    ABC_new == "" ~ NA, 
                                    TRUE ~ NA), 
                Braak_cat = case_when(Braak_Stage_new %in% c("Stage 0", "Stage 1", "Stage 2", "Stage 3") ~ "Lower stage", 
                                      Braak_Stage_new %in% c("Stage 4", "Stage 5", "Stage 6") ~ "Higher stage", 
                                      Braak_Stage_new == "" ~ NA, 
                                      TRUE ~ NA), 
                CERAD_cat = case_when(CERAD_new %in% c("Frequent", "Moderate") ~ "Frequent", 
                                      CERAD_new %in% c("No", "Sparse") ~ "Sparse", 
                                      TRUE ~ NA))

# 3- extract pairs of non-mixture --------------
pairs %>% 
  dplyr::filter(!str_detect(Exposure, "mixture")) %>% 
  dplyr::select(mz, time, Components, Mode, Exposure, Window, Pvalue_a, P_fdr) %>% 
  left_join(metabolites, by = c("mz", "time", "Components", "Mode")) %>%
  dplyr::filter(!is.na(metabolite_1)) %>% 
  unique() -> pair_1
# pair_1 %>% View()

###test block###
i = 1
strat = "ApoE_tri"
pair_1$Window[1]
################

# 3.1 Metabolomics-wide Association Study (MWAS) for individual pollutant
Ex.coef = function(i, strat){

  # assign data
  if(pair_1$Mode[i] == "Positive" & pair_1$Window[i] == "10y"){
    data = df.hilic %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_", pair_1$Exposure[i], "_10y")]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Negative" & pair_1$Window[i] == "10y"){
    data = df.c18 %>% 
      dplyr::filter(!is.na(df.c18[[str_c("SEDAC_", pair_1$Exposure[i], "_10y")]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Positive" & pair_1$Window[i] != "10y"){
    data = df.hilic %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_", pair_1$Exposure[i], "_", pair_1$Window[i])]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Negative" & pair_1$Window[i] != "10y"){
    data = df.c18 %>% 
      dplyr::filter(!is.na(df.c18[[str_c("SEDAC_", pair_1$Exposure[i], "_", pair_1$Window[i])]]) & !is.na(V1))
  }

  print(glue("For TRAP [ {pair_1$Exposure[i]} ], a total of [ {nrow(data)} ] obs were obtained for downstream analysis!")) 

  # MWAS equations
  switch(strat,
  "sex" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # stratification analyses model where Age_MY was stratified and removed from the model
  },
  "Age_cat65" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # stratification analyses model where Age_MY was stratified and removed from the model
  },
  "Age_cat75" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # stratification analyses model where Age_MY was stratified and removed from the model
  },
  "Age_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # stratification analyses model where Age_MY was stratified and removed from the model
  },
  "ApoE_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK") # stratification analyses model where ApoE_tri was stratified and removed from the model
  }, 
  "ApoE_tri" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK") # stratification analyses model where ApoE_tri was stratified and removed from the model
  }, 
  "ABC_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # stratification analyses model (same one)
    },
  "Braak_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # stratification analyses model (same one)
  },
  "CERAD_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # stratification analyses model (same one)
  }  
  )

  # split datasets into subdata according to strat
  if(strat == "Age_cat"){
    data.1 = data %>% dplyr::filter(Age_cat == "<70 yr")
    data.2 = data %>% dplyr::filter(Age_cat == ">=70 yr")
    strat.1 = "<70 yr"
    strat.2 = ">=70 yr"
  } else if(strat == "sex"){
    data.1 = data %>% dplyr::filter(sex == "Male")
    data.2 = data %>% dplyr::filter(sex == "Female")
    strat.1 = "Male"
    strat.2 = "Female"
  } else if(strat == "Age_cat65"){
    data.1 = data %>% dplyr::filter(Age_cat65 == "<65 yr")
    data.2 = data %>% dplyr::filter(Age_cat65 == ">=65 yr")
    strat.1 = "<65 yr"
    strat.2 = ">=65 yr"
  } else if(strat == "Age_cat75"){
    data.1 = data %>% dplyr::filter(Age_cat75 == "<75 yr")
    data.2 = data %>% dplyr::filter(Age_cat75 == ">=75 yr")
    strat.1 = "<75 yr"
    strat.2 = ">=75 yr"
  } else if(strat == "ApoE_cat"){
    data.1 = data %>% dplyr::filter(ApoE_cat == "E4 absence") 
    data.2 = data %>% dplyr::filter(ApoE_cat == "E4 existance")
    strat.1 = "E4 absence"
    strat.2 = "E4 existance"
  } else if(strat == "ApoE_tri"){
    data.1 = data %>% dplyr::filter(ApoE_tri == "E4 absence") 
    data.2 = data %>% dplyr::filter(ApoE_tri == "E4 heterozygous")
    data.3 = data %>% dplyr::filter(ApoE_tri == "E4 homozygous")
    strat.1 = "E4 absence"
    strat.2 = "E4 heterozygous"
    strat.3 = "E4 homozygous"
  } else if(strat == "ABC_cat"){
    data.1 = data %>% dplyr::filter(ABC_cat == "Low") 
    data.2 = data %>% dplyr::filter(ABC_cat == "High")
    strat.1 = "Low"
    strat.2 = "High"
  } else if(strat == "Braak_cat"){
    data.1 = data %>% dplyr::filter(Braak_cat == "Lower stage") 
    data.2 = data %>% dplyr::filter(Braak_cat == "Higher stage")
    strat.1 = "Lower stage"
    strat.2 = "Higher stage"
  } else if(strat == "CERAD_cat"){
    data.1 = data %>% dplyr::filter(CERAD_cat == "Sparse") 
    data.2 = data %>% dplyr::filter(CERAD_cat == "Frequent")
    strat.1 = "Sparse"
    strat.2 = "Frequent"
  } 

  # fit model for main dataset
  IQR = quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.75) - quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.25)
  
  fit0 = lm(as.formula(coding0), data = data)
  tmp <- summary(fit0)$coefficients
  if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
    # in some cases, model may not fit (due to low variance of specific metabolites)
    # so we need to assign tmp1 as follow:
    tmp0 <- cbind(n_obs = nrow(data),
                  aCoef = NA,
                  LCI_a = NA,
                  UCI_a = NA, 
                  t_a = NA,
                  Pvalue_a = NA)
    print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
  }else{      
    tmp0 <- cbind(n_obs = nrow(data),
                  aCoef = coef(fit0)[[2]]*IQR,
                  LCI_a = confint.default(fit0)[2,1]*IQR,
                  UCI_a = confint.default(fit0)[2,2]*IQR, 
                  t_a = tmp[2,3],
                  Pvalue_a = tmp[2,4])
  }      

  # fit model for data.1 -----
  IQR = quantile(as.matrix(data.1[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.75) - quantile(as.matrix(data.1[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.25)

  fit1 = lm(as.formula(coding1), data = data.1)
  tmp <- summary(fit1)$coefficients
  if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
    # in some cases, model may not fit (due to low variance of specific metabolites)
    # so we need to assign tmp1 as follow:
    tmp1 <- cbind(n_obs = nrow(data.1),
                  aCoef = NA,
                  LCI_a = NA,
                  UCI_a = NA, 
                  t_a = NA,
                  Pvalue_a = NA)
    print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
  }else{      
    tmp1 <- cbind(n_obs = nrow(data.1),
                  aCoef = coef(fit1)[[2]]*IQR,
                  LCI_a = confint.default(fit1)[2,1]*IQR,
                  UCI_a = confint.default(fit1)[2,2]*IQR, 
                  t_a = tmp[2,3],
                  Pvalue_a = tmp[2,4])
  }

  # fit model for data.2 -----
  IQR = quantile(as.matrix(data.2[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.75) - quantile(as.matrix(data.2[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.25)

  fit2 = lm(as.formula(coding1), data = data.2)
  tmp <- summary(fit2)$coefficients
  if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
    # in some cases, model may not fit (due to low variance of specific metabolites)
    # so we need to assign tmp1 as follow:
    tmp2 <- cbind(n_obs = nrow(data.2),
                  aCoef = NA,
                  LCI_a = NA,
                  UCI_a = NA, 
                  t_a = NA,
                  Pvalue_a = NA)
    print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
  }else{      
    tmp2 <- cbind(n_obs = nrow(data.2),
                  aCoef = coef(fit2)[[2]]*IQR,
                  LCI_a = confint.default(fit2)[2,1]*IQR,
                  UCI_a = confint.default(fit2)[2,2]*IQR, 
                  t_a = tmp[2,3],
                  Pvalue_a = tmp[2,4])
  }  

  if(strat == "ApoE_tri"){
    # fit model for data.3 -----
    IQR = quantile(as.matrix(data.3[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.75) - quantile(as.matrix(data.3[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.25)

    fit3 = lm(as.formula(coding1), data = data.3)
    tmp <- summary(fit3)$coefficients
    if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
      # in some cases, model may not fit (due to low variance of specific metabolites)
      # so we need to assign tmp1 as follow:
      tmp3 <- cbind(n_obs = nrow(data.3),
                    aCoef = NA,
                    LCI_a = NA,
                    UCI_a = NA, 
                    t_a = NA,
                    Pvalue_a = NA)
      print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
    }else{      
      tmp3 <- cbind(n_obs = nrow(data.3),
                    aCoef = coef(fit3)[[2]]*IQR,
                    LCI_a = confint.default(fit3)[2,1]*IQR,
                    UCI_a = confint.default(fit3)[2,2]*IQR, 
                    t_a = tmp[2,3],
                    Pvalue_a = tmp[2,4])
    }  
  }

  if(strat == "ApoE_tri"){
      # combine coefficients
    rbind(tmp0 %>% 
            as_tibble() %>% 
            dplyr::mutate(Strat_grp = strat, Strat = "All", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")),
          tmp1 %>% 
            as_tibble() %>% 
            dplyr::mutate(Strat_grp = strat, Strat = strat.1, Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")), 
          tmp2 %>% 
            as_tibble() %>% 
            dplyr::mutate(Strat_grp = strat, Strat = strat.2, Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")), 
          tmp3 %>% 
            as_tibble() %>% 
            dplyr::mutate(Strat_grp = strat, Strat = strat.3, Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}"))  ) %>% 
        dplyr::mutate(Components = glue("{pair_1$Components[i]}"), 
                      Mode = pair_1$Mode[i], 
                      Metabolite = pair_1$metabolite_1[i]) %>% 
        dplyr::select(Strat_grp, Strat, Exposure, Components, Mode, Metabolite, everything()) -> tmp.mg
  } else {
    # combine coefficients
    rbind(tmp0 %>% 
          as_tibble() %>% 
          dplyr::mutate(Strat_grp = strat, Strat = "All", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")),
        tmp1 %>% 
          as_tibble() %>% 
          dplyr::mutate(Strat_grp = strat, Strat = strat.1, Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")), 
        tmp2 %>% 
          as_tibble() %>% 
          dplyr::mutate(Strat_grp = strat, Strat = strat.2, Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}"))) %>% 
      dplyr::mutate(Components = glue("{pair_1$Components[i]}"), 
                    Mode = pair_1$Mode[i], 
                    Metabolite = pair_1$metabolite_1[i]) %>% 
      dplyr::select(Strat_grp, Strat, Exposure, Components, Mode, Metabolite, everything()) -> tmp.mg
  }
  return(tmp.mg)

}

# 3.2 execute individual pollutant (pair_1), excluding mixture pairs
# sex
tictoc::tic()
strat.res_sex = purrr::map2_dfr(1:nrow(pair_1),
                rep("sex", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
strat.res_sex %>% 
  fwrite(str_c(outpath, "/DataAnalysis_4_Stratification", 
              "/SEDAC_individual_pollutant_stratified_by_sex.csv"))


# 3- extract pairs of non-mixture --------------
pairs %>% 
  dplyr::filter(str_detect(Exposure, "mixture")) %>% 
  dplyr::select(mz, time, Components, Mode, Exposure, Window, Pvalue_a, P_fdr) %>% 
  left_join(metabolites, by = c("mz", "time", "Components", "Mode")) %>%
  dplyr::filter(!is.na(metabolite_1)) %>% 
  unique() -> pair_2

################
# since pair_2 is empty, indicating no Level-1 metabolites have been identified associated with mixture exposure (FDR q<0.2)
################

################
############# interaction analysis -------------------------
################
rm(list=(ls()))

# 1- set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# load data
load(str_c(outpath, "/!Data_for_Data_Analysis_1.rdata"))
ls()

# read metabolite list to be used
metabolites = fread(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_fdr0.2_final.csv")) %>% 
              dplyr::select(mz, time, Components, Mode, metabolite_1)

# read pairs need to be modeled for stratification analyses
pairs = fread(str_c(outpath, "/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
                                 "/Annotated_metabolites_XtoM_SEDAC_thresholds_mz10ppm_fdr0.2.csv"))


# 2- features used for stratification --------------------
df.c18 %<>% 
  dplyr::mutate(
                Age_cat65 = case_when(Age_MY < 65 ~ "<65 yr", 
                                    Age_MY >= 65 ~ ">=65 yr", 
                                    TRUE ~ NA),     
                Age_cat75 = case_when(Age_MY < 75 ~ "<75 yr", 
                                    Age_MY >= 75 ~ ">=75 yr", 
                                    TRUE ~ NA),     
                Age_cat = case_when(Age_MY < 70 ~ "<70 yr", 
                                    Age_MY >= 70 ~ ">=70 yr", 
                                    TRUE ~ NA), 
                ApoE_cat = case_when(ApoE_tri  %in% c("E4 homozygous", "E4 heterozygous") ~ "E4 existance", 
                                     ApoE_tri == "E4 absence" ~ "E4 absence", 
                                     TRUE ~ NA), 
                ABC_cat = case_when(ABC_new %in% c("High", "Intermediate") ~ "High", 
                                    ABC_new %in% c("Low", "Not") ~ "Low",
                                    ABC_new == "" ~ NA, 
                                    TRUE ~ NA), 
                Braak_cat = case_when(Braak_Stage_new %in% c("Stage 0", "Stage 1", "Stage 2", "Stage 3") ~ "Lower stage", 
                                      Braak_Stage_new %in% c("Stage 4", "Stage 5", "Stage 6") ~ "Higher stage", 
                                      Braak_Stage_new == "" ~ NA, 
                                      TRUE ~ NA), 
                CERAD_cat = case_when(CERAD_new %in% c("Frequent", "Moderate") ~ "Frequent", 
                                      CERAD_new %in% c("No", "Sparse") ~ "Sparse", 
                                      TRUE ~ NA))
df.hilic %<>% 
  dplyr::mutate(
                Age_cat65 = case_when(Age_MY < 65 ~ "<65 yr", 
                                    Age_MY >= 65 ~ ">=65 yr", 
                                    TRUE ~ NA),     
                Age_cat75 = case_when(Age_MY < 75 ~ "<75 yr", 
                                    Age_MY >= 75 ~ ">=75 yr", 
                                    TRUE ~ NA), 
                Age_cat = case_when(Age_MY < 70 ~ "<70 yr", 
                                    Age_MY >= 70 ~ ">=70 yr", 
                                    TRUE ~ NA), 
                ApoE_cat = case_when(ApoE_tri  %in% c("E4 homozygous", "E4 heterozygous") ~ "E4 existance", 
                                     ApoE_tri == "E4 absence" ~ "E4 absence", 
                                     TRUE ~ NA), 
                ABC_cat = case_when(ABC_new %in% c("High", "Intermediate") ~ "High", 
                                    ABC_new %in% c("Low", "Not") ~ "Low",
                                    ABC_new == "" ~ NA, 
                                    TRUE ~ NA), 
                Braak_cat = case_when(Braak_Stage_new %in% c("Stage 0", "Stage 1", "Stage 2", "Stage 3") ~ "Lower stage", 
                                      Braak_Stage_new %in% c("Stage 4", "Stage 5", "Stage 6") ~ "Higher stage", 
                                      Braak_Stage_new == "" ~ NA, 
                                      TRUE ~ NA), 
                CERAD_cat = case_when(CERAD_new %in% c("Frequent", "Moderate") ~ "Frequent", 
                                      CERAD_new %in% c("No", "Sparse") ~ "Sparse", 
                                      TRUE ~ NA))

# 3- extract pairs of non-mixture --------------
pairs %>% 
  dplyr::filter(!str_detect(Exposure, "mixture")) %>% 
  dplyr::select(mz, time, Components, Mode, Exposure, Window, Pvalue_a, P_fdr) %>% 
  left_join(metabolites, by = c("mz", "time", "Components", "Mode")) %>%
  dplyr::filter(!is.na(metabolite_1)) %>% 
  unique() -> pair_1


# 3.1 Metabolomics-wide Association Study (MWAS) for individual pollutant, accounting for interaction terms
Ex.inter = function(i, strat){

  # assign data
  if(pair_1$Mode[i] == "Positive" & pair_1$Window[i] == "10y"){
    data = df.hilic %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_", pair_1$Exposure[i], "_10y")]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Negative" & pair_1$Window[i] == "10y"){
    data = df.c18 %>% 
      dplyr::filter(!is.na(df.c18[[str_c("SEDAC_", pair_1$Exposure[i], "_10y")]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Positive" & pair_1$Window[i] != "10y"){
    data = df.hilic %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_", pair_1$Exposure[i], "_", pair_1$Window[i])]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Negative" & pair_1$Window[i] != "10y"){
    data = df.c18 %>% 
      dplyr::filter(!is.na(df.c18[[str_c("SEDAC_", pair_1$Exposure[i], "_", pair_1$Window[i])]]) & !is.na(V1))
  }

  print(glue("For TRAP [ {pair_1$Exposure[i]} ], a total of [ {nrow(data)} ] obs were obtained for downstream analysis!")) 

  # MWAS equations
  switch(strat, 
  "sex" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*sex + Age_MY + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analysis model where Age_cat was interacted with pollutants
  },
  "Age_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*Age_cat + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analysis model where Age_cat was interacted with pollutants
  },
  "Age_cat65" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*Age_cat65 + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analysis model where Age_cat was interacted with pollutants
  },
  "Age_cat75" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*Age_cat75 + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analysis model where Age_cat was interacted with pollutants
  },
  "ApoE_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*ApoE_cat + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK") # interaction analyses model where ApoE_tri was stratified and removed from the model
  }, 
  "ApoE_tri" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*ApoE_tri + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK") # interaction analyses model where ApoE_tri was stratified and removed from the model
  }, 
  "ABC_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*ABC_cat + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analyses model (same one)
    },
  "Braak_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*Braak_cat + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analyses model (same one)
  },
  "CERAD_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*CERAD_cat + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analyses model (same one)
  }  
  )

  # fit model for main dataset
  IQR = quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.75) - quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.25)
  
  fit0 = lm(as.formula(coding0), data = data)
  tmp <- summary(fit0)$coefficients
  if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
    # in some cases, model may not fit (due to low variance of specific metabolites)
    # so we need to assign tmp1 as follow:
    tmp0 <- cbind(n_obs = nrow(data),
                  aCoef = NA,
                  LCI_a = NA,
                  UCI_a = NA, 
                  t_a = NA,
                  Pvalue_a = NA)
    print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
  }else{      
    tmp0 <- cbind(n_obs = nrow(data),
                  aCoef = coef(fit0)[[2]]*IQR,
                  LCI_a = confint.default(fit0)[2,1]*IQR,
                  UCI_a = confint.default(fit0)[2,2]*IQR, 
                  t_a = tmp[2,3],
                  Pvalue_a = tmp[2,4])
  }      

  # fit interaction model -----
  IQR = quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.75) - quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.25)

  fit1 = lm(as.formula(coding1), data = data)
  tmp <- summary(fit1)$coefficients

  if(strat == "ApoE_tri"){
    if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
      # in some cases, model may not fit (due to low variance of specific metabolites)
      # so we need to assign tmp1 as follow:
      tmp1.1 <- cbind(n_obs = nrow(data),
                    aCoef = NA,
                    LCI_a = NA,
                    UCI_a = NA, 
                    t_a = NA,
                    Pvalue_a = NA) # interaction model--main effect
      tmp1.2 <- cbind(n_obs = nrow(data),
                    aCoef = NA,
                    LCI_a = NA,
                    UCI_a = NA, 
                    t_a = NA,
                    Pvalue_a = NA) # interaction model--main effect
      tmp2.1 <- cbind(n_obs = nrow(data),
                    aCoef = NA,
                    LCI_a = NA,
                    UCI_a = NA, 
                    t_a = NA,
                    Pvalue_a = NA) # interaction model--interaction effect
      tmp2.2 <- cbind(n_obs = nrow(data),
                    aCoef = NA,
                    LCI_a = NA,
                    UCI_a = NA, 
                    t_a = NA,
                    Pvalue_a = NA) # interaction model--interaction effect
      print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
    }else{      
      tmp1.1 <- cbind(n_obs = nrow(data),
                      aCoef = coef(fit1)[[2]]*IQR,
                      LCI_a = confint.default(fit1)[2,1]*IQR,
                      UCI_a = confint.default(fit1)[2,2]*IQR, 
                      t_a = tmp[2,3],
                      Pvalue_a = tmp[2,4]) # interaction model--main effect
      tmp1.2 <- cbind(n_obs = nrow(data),
                      aCoef = coef(fit1)[[3]]*IQR,
                      LCI_a = confint.default(fit1)[3,1]*IQR,
                      UCI_a = confint.default(fit1)[3,2]*IQR, 
                      t_a = tmp[3,3],
                      Pvalue_a = tmp[3,4]) # interaction model--main effect
      tmp2.1 <- cbind(n_obs = nrow(data),
                      aCoef = coef(fit1)[[length(coef(fit1))-1]]*IQR,
                      LCI_a = confint.default(fit1)[nrow(confint.default(fit1))-1,1]*IQR,
                      UCI_a = confint.default(fit1)[nrow(confint.default(fit1))-1,2]*IQR, 
                      t_a = tmp[nrow(tmp)-1,3],
                      Pvalue_a = tmp[nrow(tmp)-1,4]) # interaction model--interaction effect
      tmp2.2 <- cbind(n_obs = nrow(data),
                      aCoef = coef(fit1)[[length(coef(fit1))]]*IQR,
                      LCI_a = confint.default(fit1)[nrow(confint.default(fit1)),1]*IQR,
                      UCI_a = confint.default(fit1)[nrow(confint.default(fit1)),2]*IQR, 
                      t_a = tmp[nrow(tmp),3],
                      Pvalue_a = tmp[nrow(tmp),4]) # interaction model--interaction effect
    }
    # combine coefficients
    rbind(tmp0 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Main analysis", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")),
          tmp1.1 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Interaction model-main effect.1", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")), 
          tmp1.2 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Interaction model-main effect.2", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")),           
          tmp2.1 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Interaction model-interaction effect.1", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")),          
          tmp2.2 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Interaction model-interaction effect.2", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}"))
          ) %>% 
        dplyr::mutate(Components = glue("{pair_1$Components[i]}"), 
                      Mode = pair_1$Mode[i], 
                      Metabolite = pair_1$metabolite_1[i]) %>% 
        dplyr::select(Interaction_term, Group, Exposure, Components, Mode, Metabolite, everything()) -> tmp.mg
  }else{
    if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
      # in some cases, model may not fit (due to low variance of specific metabolites)
      # so we need to assign tmp1 as follow:
      tmp1 <- cbind(n_obs = nrow(data),
                    aCoef = NA,
                    LCI_a = NA,
                    UCI_a = NA, 
                    t_a = NA,
                    Pvalue_a = NA) # interaction model--main effect
      tmp2 <- cbind(n_obs = nrow(data),
                    aCoef = NA,
                    LCI_a = NA,
                    UCI_a = NA, 
                    t_a = NA,
                    Pvalue_a = NA) # interaction model--interaction effect
      print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
    }else{      
      tmp1 <- cbind(n_obs = nrow(data),
                    aCoef = coef(fit1)[[2]]*IQR,
                    LCI_a = confint.default(fit1)[2,1]*IQR,
                    UCI_a = confint.default(fit1)[2,2]*IQR, 
                    t_a = tmp[2,3],
                    Pvalue_a = tmp[2,4]) # interaction model--main effect
      tmp2 <- cbind(n_obs = nrow(data),
                    aCoef = coef(fit1)[[length(coef(fit1))]]*IQR,
                    LCI_a = confint.default(fit1)[nrow(confint.default(fit1)),1]*IQR,
                    UCI_a = confint.default(fit1)[nrow(confint.default(fit1)),2]*IQR, 
                    t_a = tmp[nrow(tmp),3],
                    Pvalue_a = tmp[nrow(tmp),4]) # interaction model--interaction effect
    }
    # combine coefficients
    rbind(tmp0 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Main analysis", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")),
          tmp1 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Interaction model-main effect", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")), 
          tmp2 %>% 
            as_tibble() %>% 
            dplyr::mutate(Interaction_term = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}*{strat}"), Group = "Interaction model-interaction effect", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}"))
          ) %>% 
        dplyr::mutate(Components = glue("{pair_1$Components[i]}"), 
                      Mode = pair_1$Mode[i], 
                      Metabolite = pair_1$metabolite_1[i]) %>% 
        dplyr::select(Interaction_term, Group, Exposure, Components, Mode, Metabolite, everything()) -> tmp.mg
  }
  return(tmp.mg)
}

# 3.2 execute individual pollutant (pair_1), excluding mixture pairs 
# sex
tictoc::tic()
strat.res_sex = purrr::map2_dfr(1:nrow(pair_1),
                rep("sex", nrow(pair_1)),
                Ex.inter)
tictoc::toc()
strat.res_sex %>% 
  fwrite(str_c(outpath, "/DataAnalysis_5_Interaction", 
              "/SEDAC_individual_pollutant_interaction_with_sex.csv"))

# 3- extract pairs of non-mixture --------------
pairs %>% 
  dplyr::filter(str_detect(Exposure, "mixture")) %>% 
  dplyr::select(mz, time, Components, Mode, Exposure, Window, Pvalue_a, P_fdr) %>% 
  left_join(metabolites, by = c("mz", "time", "Components", "Mode")) %>%
  dplyr::filter(!is.na(metabolite_1)) %>% 
  unique() -> pair_2

################
# since pair_2 is empty, indicating no Level-1 metabolites have been identified associated with mixture exposure (FDR q<0.2)
################



# viz ---------------------
library(patchwork)
library(RColorBrewer)

## sex ----
fread(str_c(outpath, "/DataAnalysis_4_Stratification",
            "/SEDAC_individual_pollutant_stratified_by_sex.csv")) %>% 
  dplyr::select(Strat, Exposure, Components, Mode, aCoef, LCI_a, UCI_a, Pvalue_a) %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"), 
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  arrange(Window, Exposure, Components, Mode) -> df2

# read metabolites information
metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))
metabolite %<>% 
  dplyr::select(Mode, Components, "Brief name") %>% 
  unique() %>% 
  dplyr::rename(Metabolite = "Brief name")
# combine
df2 %<>% 
  dplyr::left_join(metabolite, by = c("Components", "Mode"))

# pairs to be used
df2$Exposure %>% unique()
df2$Window %>% unique()
df2$Metabolite %>% unique()
# reset levels arranged by metabolites
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Glucosamine",
              "Vitamin D2")
# re-assign levels
df2 %>% 
  dplyr::mutate(Pair = str_c(Window, Exposure, Metabolite, sep = "_"), 
                Strat_grp = factor(Strat, levels = c("Male", "Female", "All")),
                Exposure = factor(Exposure, levels = c("O3", "SO4", "Silicon", "Lead", "Vanadium")), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y")), 
                Metabolite = factor(Metabolite, levels = meta.lvl)) %>% 
  arrange(Metabolite, Window, Exposure, Strat_grp) -> plt.dt.2
plt.dt.2$Pair %>% unique() -> pair.lvl
plt.dt.2 %<>% 
  dplyr::mutate(Pair = factor(Pair, levels = pair.lvl))

# import interaction analysis results
fread(str_c(outpath, "/DataAnalysis_5_Interaction",
            "/SEDAC_individual_pollutant_interaction_with_sex.csv")) %>% 
  dplyr::filter(Group == "Interaction model-interaction effect") %>% 
  dplyr::filter(Pvalue_a <0.05) %>% 
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"),
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::select(Exposure, Window, Components, Mode) %>% 
  dplyr::mutate(flag="*", Strat_grp="All") -> inter
# combine interaction info
plt.dt.2 %<>% 
  left_join(inter, by = c("Exposure", "Window", "Components", "Mode", "Strat_grp"))

plt.dt.2 %<>% 
  dplyr::mutate(Strat_grp = factor(Strat_grp, levels = c("Male", "Female", "All")))

# sub-figures 
# colors
meta.color = c(
  "#4A6FA5", "#C97A7A", "#6E7F5F", "#9A85C9", "#D9A066",
  "#8AB4D0", "#B85C57", "#9BAF7F", "#7765A8", "#E8BB7F",
  "#6B8FBF", "#8C3F3F", "#556B2F", "#C2AEDB", "#F0D9A8",
  "#A49E8D", "#69645C"
)

p1 = plt.dt.2 %>% # time windows
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Window)) + 
  scale_fill_manual(values = c("#c0bdd1", "#8881BC", "#4850A0", "#103B95")) + 
  theme_minimal() +
  theme(#axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p2 = plt.dt.2 %>% # exposures
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Exposure)) +
  scale_fill_manual(values = c("#e58080", "#6fcab4", "#7f5b9c", "#a5af4c", "#a1a19f")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p3 = plt.dt.2 %>% 
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Metabolite)) +
  scale_fill_manual(values = meta.color) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p4 = plt.dt.2 %>%
  ggplot(aes(x = aCoef, y = Pair)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "#525050") +
  geom_point(aes(color = Strat_grp),
             size = 3,
             position = position_dodge(0.7)) + 
  geom_errorbar(aes(xmin = LCI_a, xmax = UCI_a, color = Strat_grp),
                width = 0.4,
                position = position_dodge(0.7)) + 
  scale_color_manual(values = c("#909e3d", "#4d5bac", "#46463e")) +
  # scale_x_continuous(limits = c(-6.5, 6), breaks = c(-5, -2.5, 0, 2.5, 5)) +
  geom_text(data = plt.dt.2,  # remember change datasets here --------------------
             aes(UCI_a+0.5, Pair, group = Strat_grp, label = flag), 
             size = 8, position = position_dodge(0.5), color = "red") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position = "right") + 
  facet_grid(.~Mode, scales = "free_x") + 
  guides(color = guide_legend(override.aes = list(size = 12)))

# merge
(p1|p2|p3|p4) + plot_layout(widths = c(0.5,0.5,0.5,7), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right",
        legend.key.size = unit(0.9, "cm"),
        axis.text = element_text(size = 18))
# save
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig4_Sex_stratum", ".png"), 
       height = 20, width = 20, dpi = 600)




pair_1$Exposure %>% View()

############################
# extract SD value for each metabolite under specified sample size
############################
i=1
# all popu
df.c18 %>% 
      dplyr::filter(!is.na(df.c18[[str_c("SEDAC_", pair_1$Exposure[11], "_10y")]]) & !is.na(V1)) %>% 
  .$V12152 %>% sd(na.rm = T)

