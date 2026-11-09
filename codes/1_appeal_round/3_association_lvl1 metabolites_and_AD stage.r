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
library(MASS) # for ordinal regression
library(brant)
library(nnet)
# 1- set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# load data
load(str_c(outpath, "/!Data_for_Data_Analysis_1.rdata"))
ls()

# ###############################
# # export air pollutant data for ArcGIS 
# fread(file = str_c(inpath, "/Demographics/Brain_bank_Marla_clean_432obs.csv")) -> demo
# # View(demo)
# demo %>% 
#   dplyr::select(Case_Number, POINT_X, POINT_Y) %>% 
#   set_names(c("Case_Number", "lon", "lat"))  -> demo.ids

# df.c18 %>% 
#   dplyr::filter(!is.na(V1) & !is.na(SEDAC_PM25NH4_1y)) %>%
#   left_join(demo.ids, by = "Case_Number") %>% 
#   dplyr::select(lon, lat, SEDAC_PM25_1y, SEDAC_NO2_1y, SEDAC_O3_1y) -> out

# out %>% 
#   writexl::write_xlsx(str_c(outpath, "/Air_pollutant_data_for_ArcGIS.xlsx"))
# ###############################


# read metabolite list to be used
metabolites <- fread(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_fdr0.2_final.csv")) %>% 
              dplyr::select(mz, time, Components, Mode, metabolite_1)

# read pairs need to be modeled for ordinal regression on AD markers
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

# make AD markers as ordered factoral variables 

df.c18$Braak_Stage_new %>% table(useNA = "always")
df.c18$CERAD_new %>% table(useNA = "always")
df.c18$ABC_new %>% table(useNA = "always")

df.c18 %<>% 
  dplyr::mutate(Braak = case_when(Braak_Stage_new == "" ~ NA, 
                                  TRUE ~ Braak_Stage_new), 
                CERAD = case_when(CERAD_new == "" ~ NA, 
                                  TRUE ~ CERAD_new), 
                ABC = case_when(ABC_new == "" ~ NA, 
                                TRUE ~ ABC_new)) %>%
  dplyr::mutate(Braak = factor(Braak, levels = c("Stage 0", "Stage 1", "Stage 2", "Stage 3", "Stage 4", "Stage 5", "Stage 6"), ordered = TRUE),
                CERAD = factor(CERAD, levels = c("No", "Sparse", "Moderate", "Frequent"), ordered = TRUE),
                ABC = factor(ABC, levels = c("Not", "Low", "Intermediate", "High"), ordered = TRUE))
df.hilic %<>% 
  dplyr::mutate(Braak = case_when(Braak_Stage_new == "" ~ NA, 
                                  TRUE ~ Braak_Stage_new), 
                CERAD = case_when(CERAD_new == "" ~ NA, 
                                  TRUE ~ CERAD_new), 
                ABC = case_when(ABC_new == "" ~ NA, 
                                TRUE ~ ABC_new)) %>%
  dplyr::mutate(Braak = factor(Braak, levels = c("Stage 0", "Stage 1", "Stage 2", "Stage 3", "Stage 4", "Stage 5", "Stage 6"), ordered = TRUE),
                CERAD = factor(CERAD, levels = c("No", "Sparse", "Moderate", "Frequent"), ordered = TRUE),
                ABC = factor(ABC, levels = c("Not", "Low", "Intermediate", "High"), ordered = TRUE))

# 3- extract pairs of non-mixture --------------
pairs %>% 
  dplyr::filter(!str_detect(Exposure, "mixture")) %>% 
  dplyr::select(mz, time, Components, Mode, Exposure, Window, Pvalue_a, P_fdr) %>% 
  left_join(metabolites, by = c("mz", "time", "Components", "Mode")) %>%
  dplyr::filter(!is.na(metabolite_1)) %>% 
  unique() -> pair_1

# only keep unique metabolites since we are working on investigating associations between these metabolites and AD markers
pair_1 %<>% 
  dplyr::select(mz, time, Components, Mode, metabolite_1) %>% 
  unique()


###test block###
i = 1
strat = "Braak"
################


# 3.1 Metabolomics-wide Association Study (MWAS) for individual pollutant
Ex.coef = function(i, strat){

  # assign data
  if(pair_1$Mode[i] == "Positive"){
    data = df.hilic %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Negative"){
    data = df.c18 %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1))
  }

  print(glue("For [ {pair_1$metabolite_1[i]} ], a total of [ {nrow(data)} ] obs were obtained for downstream analysis!")) 

  # MWAS equations
  switch(strat,
  "Braak" = {
    coding0 = glue("Braak ~ {pair_1$Components[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
  },
  "ABC" = {
    coding0 = glue("ABC ~ {pair_1$Components[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
  },
  "CERAD" = {
    coding0 = glue("CERAD ~ {pair_1$Components[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
  }
  )

  # fit model for main dataset
  fit0 = polr(as.formula(coding0), data = data, Hess = TRUE)
  fit1 = multinom(as.formula(coding0), data = data)

  # does relaxing propotional odds improve fit? 
  LR <- deviance(fit0) - deviance(fit1)
  df <- fit1$edf - fit0$edf
  po_p <-pchisq(LR, df, lower.tail = FALSE)

  tmp <- summary(fit0)$coefficients
  
  if(is.null(tmp) | !(glue("{pair_1$Components[i]}") %in% row.names(tmp))){
    # in some cases, model may not fit (due to low variance of specific metabolites)
    # so we need to assign tmp1 as follow:
    tmp0 <- cbind(n_obs = nrow(data),
                  aCoef = NA,
                  LCI_a = NA,
                  UCI_a = NA, 
                  t_a = NA,
                  Pvalue_a = NA, 
                  PO_p = NA)
    print(str_c("Warning!! metabolite [", glue("{pair_1$Components[i]}"), "] could not fit adjusted models for ", glue("SEDAC_{pair_1$Exposure[i]}_{pair_1$Window[i]}"), " due to low variates !"))
  }else{      
    tmp0 <- cbind(n_obs = nrow(data),
                  aCoef = coef(fit0)[[1]],
                  LCI_a = confint.default(fit0)[1,1],
                  UCI_a = confint.default(fit0)[1,2], 
                  t_a = tmp[1,3],
                  Pvalue_a = 2*(1-pnorm(abs(tmp[1,3]))),
                  PO_p = po_p)
  }      
  tmp0 %>% 
    as_tibble() %>% 
    dplyr::mutate(Components = glue("{pair_1$Components[i]}"), 
                  Mode = pair_1$Mode[i], 
                  Metabolite = pair_1$metabolite_1[i]) %>% 
      dplyr::select(Components, Mode, Metabolite, everything()) -> res
  return(res)

}

# 3.2 execute individual pollutant (pair_1), excluding mixture pairs
# Braak
tictoc::tic()
res_Braak = purrr::map2_dfr(1:nrow(pair_1),
                rep("Braak", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
res_Braak %>% 
  fwrite(str_c(outpath, "/DataAnalysis_7_association_with_AD", 
              "/Association_between_lvl1_metabolites_and_Braak.csv"))
res_Braak %>% View()

# ABC
tictoc::tic()
res_ABC = purrr::map2_dfr(1:nrow(pair_1),
                rep("ABC", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
res_ABC %>% 
  fwrite(str_c(outpath, "/DataAnalysis_7_association_with_AD", 
              "/Association_between_lvl1_metabolites_and_ABC.csv"))

# CERAD
tictoc::tic()
res_CERAD = purrr::map2_dfr(1:nrow(pair_1),
                rep("CERAD", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
res_CERAD %>% 
  fwrite(str_c(outpath, "/DataAnalysis_7_association_with_AD", 
              "/Association_between_lvl1_metabolites_and_CERAD.csv"))





# alternative: use binary outcomes as we did in stratification analyses -------------------
df.c18$Braak_cat %>% table(useNA = "always")
df.c18$CERAD_cat %>% table(useNA = "always")
df.c18$ABC_cat %>% table(useNA = "always")

# binary variables
df.c18 %<>% 
  dplyr::mutate(Braak_cat = factor(Braak_cat, levels = c("Lower stage", "Higher stage"), ordered = TRUE),
                CERAD_cat = factor(CERAD_cat, levels = c("Sparse", "Frequent"), ordered = TRUE),
                ABC_cat = factor(ABC_cat, levels = c("Low", "High"), ordered = TRUE))
df.hilic %<>% 
  dplyr::mutate(Braak_cat = factor(Braak_cat, levels = c("Lower stage", "Higher stage"), ordered = TRUE),
                CERAD_cat = factor(CERAD_cat, levels = c("Sparse", "Frequent"), ordered = TRUE),
                ABC_cat = factor(ABC_cat, levels = c("Low", "High"), ordered = TRUE))

# 3- extract pairs of non-mixture --------------
pair_1

###test block###
i = 1
strat = "Braak"
################

# 3.1 Metabolomics-wide Association Study (MWAS) for individual pollutant
Ex.coef = function(i, strat){

  # assign data
  if(pair_1$Mode[i] == "Positive"){
    data = df.hilic %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1))
  } else if(pair_1$Mode[i] == "Negative"){
    data = df.c18 %>% 
      dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1))
  }

  print(glue("For [ {pair_1$metabolite_1[i]} ], a total of [ {nrow(data)} ] obs were obtained for downstream analysis!")) 

  # MWAS equations
  switch(strat,
  "Braak" = {
    coding0 = glue("Braak_cat ~ {pair_1$Components[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
  },
  "ABC" = {
    coding0 = glue("ABC_cat ~ {pair_1$Components[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
  },
  "CERAD" = {
    coding0 = glue("CERAD_cat ~ {pair_1$Components[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri") # main analysis model
  }
  )

  # fit model for main dataset
  fit0 = glm(as.formula(coding0), family = binomial(link = "logit"), data = data)

  tmp <- summary(fit0)$coefficients
  
  if(is.null(tmp) | !(glue("{pair_1$Components[i]}") %in% row.names(tmp))){
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
                  aCoef = coef(fit0)[[2]],
                  LCI_a = confint.default(fit0)[2,1],
                  UCI_a = confint.default(fit0)[2,2], 
                  t_a = tmp[2,3],
                  Pvalue_a = tmp[2,4])
  }      
  tmp0 %>% 
    as_tibble() %>% 
    dplyr::mutate(Components = glue("{pair_1$Components[i]}"), 
                  Mode = pair_1$Mode[i], 
                  Metabolite = pair_1$metabolite_1[i]) %>% 
      dplyr::select(Components, Mode, Metabolite, everything()) -> res
  return(res)

}

# 3.2 execute individual pollutant (pair_1), excluding mixture pairs
# Braak
tictoc::tic()
res_Braak_binary = purrr::map2_dfr(1:nrow(pair_1),
                rep("Braak", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
res_Braak_binary %>% 
  fwrite(str_c(outpath, "/DataAnalysis_7_association_with_AD", 
              "/Association_between_lvl1_metabolites_and_Braak_binary.csv"))

# ABC
tictoc::tic()
res_ABC_binary = purrr::map2_dfr(1:nrow(pair_1),
                rep("ABC", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
res_ABC_binary %>% 
  fwrite(str_c(outpath, "/DataAnalysis_7_association_with_AD", 
              "/Association_between_lvl1_metabolites_and_ABC_binary.csv"))

# CERAD
tictoc::tic()
res_CERAD_binary = purrr::map2_dfr(1:nrow(pair_1),
                rep("CERAD", nrow(pair_1)),
                Ex.coef)
tictoc::toc()
res_CERAD_binary %>% 
  fwrite(str_c(outpath, "/DataAnalysis_7_association_with_AD", 
              "/Association_between_lvl1_metabolites_and_CERAD_binary.csv"))
