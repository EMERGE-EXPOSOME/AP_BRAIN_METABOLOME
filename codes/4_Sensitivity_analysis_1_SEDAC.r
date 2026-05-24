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


colnames(df.c18)[1:200]
df.c18$death_yr

# 2- features used for stratification --------------------
df.c18$Age_MY
df.c18$ApoE_tri %>% table(useNA = "always")
df.c18$ABC_new %>% table(useNA = "always")
df.c18$Braak_Stage_new %>% table(useNA = "always")
df.c18$CERAD_new %>% table(useNA = "always")



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
pair_1$Window[1]
cov = "death_yr" # ABC_new, Braak_Stage_new, CERAD_new
################

# 3.1 Metabolomics-wide Association Study (MWAS) for individual pollutant, accounting for interaction terms
Ex.coef = function(i, cov){

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
  switch(cov, 
  "death_yr" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + death_yr + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # sensitivity analysis by adding death_yr
  },
  "ABC_new" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + ABC_new + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # sensitivity analysis by adding ABC_new
    },
  "Braak_Stage_new" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Braak_Stage_new + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri")
  },
  "CERAD_new" = {
    coding0 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
    coding1 = glue("{pair_1$Components[i]} ~ SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]} + CERAD_new + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") 
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

  # fit sensitivity analysis model -----
  IQR = quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.75) - quantile(as.matrix(data[,glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}")]),0.25)

  fit1 = lm(as.formula(coding1), data = data)
  tmp <- summary(fit1)$coefficients
  if(is.null(tmp) | !(glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}") %in% row.names(tmp))){
    # in some cases, model may not fit (due to low variance of specific metabolites)
    # so we need to assign tmp1 as follow:
    tmp1 <- cbind(n_obs = nrow(data),
                  aCoef = NA,
                  LCI_a = NA,
                  UCI_a = NA, 
                  t_a = NA,
                  Pvalue_a = NA) # sensitivity analysis
    print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
  }else{      
    tmp1 <- cbind(n_obs = nrow(data),
                  aCoef = coef(fit1)[[2]]*IQR,
                  LCI_a = confint.default(fit1)[2,1]*IQR,
                  UCI_a = confint.default(fit1)[2,2]*IQR, 
                  t_a = tmp[2,3],
                  Pvalue_a = tmp[2,4]) # sensitivity analysis
  }

  # combine coefficients
  rbind(tmp0 %>% 
          as_tibble() %>% 
          dplyr::mutate(Group = "Main analysis", Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}")),
        tmp1 %>% 
          as_tibble() %>% 
          dplyr::mutate(Group = glue("ADd [{cov}]"), Exposure = glue("{pair_1$Exposure[i]}_{pair_1$Window[i]}"))
        ) %>% 
      dplyr::mutate(Components = glue("{pair_1$Components[i]}"), 
                    Mode = pair_1$Mode[i], 
                    Metabolite = pair_1$metabolite_1[i]) %>% 
      dplyr::select(Group, Exposure, Components, Mode, Metabolite, everything()) -> tmp3
  
  return(tmp3)

}

# 3.2 execute individual pollutant (pair_1), excluding mixture pairs 
# death_yr
tictoc::tic()
sensi.res_death_yr = purrr::map2_dfr(1:nrow(pair_1),
                     rep("death_yr", nrow(pair_1)),
                     Ex.coef)
tictoc::toc()
sensi.res_death_yr %>% 
  fwrite(str_c(outpath, "/DataAnalysis_6_Sensitivity_analysis", 
              "/SEDAC_individual_pollutant_sensitivitiy_analysis_Adding_death_yr.csv"))

# ABC_new
tictoc::tic()
sensi.res_ABC_new = purrr::map2_dfr(1:nrow(pair_1),
                     rep("ABC_new", nrow(pair_1)),
                     Ex.coef)
tictoc::toc()
sensi.res_ABC_new %>% 
  fwrite(str_c(outpath, "/DataAnalysis_6_Sensitivity_analysis", 
              "/SEDAC_individual_pollutant_sensitivitiy_analysis_Adding_ABC_new.csv"))

# Braak_Stage_new
tictoc::tic()
sensi.res_Braak_Stage_new = purrr::map2_dfr(1:nrow(pair_1),
                     rep("Braak_Stage_new", nrow(pair_1)),
                     Ex.coef)
tictoc::toc()
sensi.res_Braak_Stage_new %>% 
  fwrite(str_c(outpath, "/DataAnalysis_6_Sensitivity_analysis", 
              "/SEDAC_individual_pollutant_sensitivitiy_analysis_Adding_Braak_Stage_new.csv"))

# CERAD_new
tictoc::tic()
sensi.res_CERAD_new= purrr::map2_dfr(1:nrow(pair_1),
                     rep("CERAD_new", nrow(pair_1)),
                     Ex.coef)
tictoc::toc()
sensi.res_CERAD_new %>% 
  fwrite(str_c(outpath, "/DataAnalysis_6_Sensitivity_analysis", 
              "/SEDAC_individual_pollutant_sensitivitiy_analysis_Adding_CERAD_new.csv"))


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
