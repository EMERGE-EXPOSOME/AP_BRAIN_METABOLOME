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
metabolites = fread(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_fdr0.2_final.csv")) %>% 
              dplyr::select(mz, time, Components, Mode, metabolite_1)

# read pairs need to be modeled for stratification analyses
pairs = fread(str_c(outpath, "/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
                                 "/Annotated_metabolites_XtoM_SEDAC_thresholds_mz10ppm_fdr0.2.csv"))


# 2- features used for stratification --------------------
df.c18$Age_MY
df.c18$ApoE_tri %>% table(useNA = "always")
df.c18$ABC_new %>% table(useNA = "always")
df.c18$Braak_Stage_new %>% table(useNA = "always")
df.c18$CERAD_new %>% table(useNA = "always")

df.c18 %<>% 
  dplyr::mutate(Age_cat = case_when(Age_MY < 70 ~ "<70 yr", 
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
  dplyr::mutate(Age_cat = case_when(Age_MY < 70 ~ "<70 yr", 
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
i = 27
strat = "ApoE_tri"
pair_1$Window[1]
################

# 3.1 Metabolomics-wide Association Study (MWAS) for individual pollutant, accounting for interaction terms
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
  "Age_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*Age_cat + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # interaction analysis model where Age_cat was interacted with pollutants
  },
  "ApoE_cat" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}*ApoE_cat + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK") # interaction analyses model where ApoE_tri was stratified and removed from the model
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
# Age_cat
tictoc::tic()
strat.res_age = purrr::map2_dfr(1:nrow(pair_1),
                rep("Age_cat", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
strat.res_age %>% 
  fwrite(str_c(outpath, "/DataAnalysis_5_Interaction", 
              "/SEDAC_individual_pollutant_interaction_with_Age_cat.csv"))

# ApoE_cat
tictoc::tic()
strat.res = purrr::map2_dfr(1:nrow(pair_1),
                rep("ApoE_cat", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
strat.res %>% 
  fwrite(str_c(outpath, "/DataAnalysis_5_Interaction", 
              "/SEDAC_individual_pollutant_interaction_with_ApoE_cat.csv"))

# ABC_cat
tictoc::tic()
strat.res = purrr::map2_dfr(1:nrow(pair_1),
                rep("ABC_cat", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
strat.res %>% 
  fwrite(str_c(outpath, "/DataAnalysis_5_Interaction", 
              "/SEDAC_individual_pollutant_interaction_with_ABC_cat.csv"))

# Braak_cat
tictoc::tic()
strat.res = purrr::map2_dfr(1:nrow(pair_1),
                rep("Braak_cat", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
strat.res %>% 
  fwrite(str_c(outpath, "/DataAnalysis_5_Interaction", 
              "/SEDAC_individual_pollutant_interaction_with_Braak_cat.csv"))

# CERAD_cat
tictoc::tic()
strat.res = purrr::map2_dfr(1:nrow(pair_1),
                rep("CERAD_cat", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
strat.res %>% 
  fwrite(str_c(outpath, "/DataAnalysis_5_Interaction", 
              "/SEDAC_individual_pollutant_interaction_with_CERAD_cat.csv"))


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
