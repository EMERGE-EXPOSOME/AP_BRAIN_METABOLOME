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
               "mediation"
               )
library(MASS) # for ordinal regression
library(brant)
library(nnet)

# 1- set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# 2- load data ----------------------------
load(str_c(outpath, "/!Data_for_Data_Analysis_1.rdata"))
ls()

# read metabolite (associated with Air pollution) list to be used
meta_expo <- fread(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_fdr0.2_final.csv")) %>% 
              dplyr::select(mz, time, Components, Mode, metabolite_1) %>% 
              dplyr::rename(Metabolite = metabolite_1)

# read metabolite (associated with AD markers) list to be used
meta_outcome <- readxl::read_xlsx(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_with_ADmaker.xlsx")) %>% 
              dplyr::select(mz, rt, Metabolite, Column)
# since all features among meta_outcome came from C18, we will only extract features from C18 for features associated with AD markers
c18_info <- fread(str_c(outpath, "/log2.transformed.metabolomics_C18.csv")) %>% 
              dplyr::mutate(Components = str_c("V", row_number()))
c18_info %<>% 
  dplyr::select(mz, time, Components)
meta_outcome %<>% 
  dplyr::rename(time = rt) %>% 
  left_join(c18_info, by = c("mz", "time")) %>%
  dplyr::select(-Column) %>% 
  dplyr::mutate(Mode = "Negative")

# combine, and here we have the union set of metabolites that significantly associated with either air pollution or AD markers. 
meta_list <- rbind(meta_expo, meta_outcome)

# 3-make pairs to be used in mediation models ------------------------------------
# ad markers, 3
AD_markers <- c("Braak", "CERAD", "ABC")
# air pollutant variables, 72
AP_vars <-df.c18 %>% 
  dplyr::select(contains("SEDAC") & !contains("_15y")) %>%  colnames()

# pairs to be used in mediation analsyis 
pairs <- expand_grid(
  exposure = AP_vars,       # 72 exposures
  disease  = AD_markers,    # 3 diseases
  meta_list                 # metabolite data frame -> each row crossed in
) %>%
  mutate(triplet_id = row_number()) %>%
  relocate(triplet_id, exposure, Components, disease)


# 4- data processing before mediation --------------------
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
# since AD markers are too sparse to be moelled under mediation framework, we simplified them as binary categorical variable
df.c18 %<>% 
  dplyr::mutate(Braak_cat = factor(Braak_cat, levels = c("Lower stage", "Higher stage"), ordered = TRUE)) %>% 
  dplyr::mutate(CERAD_cat = factor(CERAD_cat, levels = c("Sparse", "Frequent"), ordered = TRUE)) %>% 
  dplyr::mutate(ABC_cat = factor(ABC_cat, levels = c("Low", "High"), ordered = TRUE))
df.hilic %<>% 
  dplyr::mutate(Braak_cat = factor(Braak_cat, levels = c("Lower stage", "Higher stage"), ordered = TRUE)) %>% 
  dplyr::mutate(CERAD_cat = factor(CERAD_cat, levels = c("Sparse", "Frequent"), ordered = TRUE)) %>% 
  dplyr::mutate(ABC_cat = factor(ABC_cat, levels = c("Low", "High"), ordered = TRUE))


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


# 5- mediation --------------------

###test block###
i = 1
################
pairs

Pairwise_1 = function(i){

  # assign data
  if(pairs$Mode[i] == "Positive"){
    data = df.hilic %>% 
      dplyr::filter(!is.na(df.hilic[[pairs$exposure[[i]]]]) & !is.na(V1))
  } else if(pairs$Mode[i] == "Negative"){
    data = df.c18 %>% 
      dplyr::filter(!is.na(df.c18[[pairs$exposure[[i]]]]) & !is.na(V1))
  }

  print(glue("For [ {pairs$Metabolite[i]} ], a total of [ {nrow(data)} ] obs were obtained for downstream analysis!")) 


  # pairwise mediation model
  coding1 = glue("{pairs$Components[i]} ~ {pairs$exposure[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri")
  fit1 = lm(as.formula(coding1), data = data)
  coding2 =  glue("{pairs$disease[i]}_cat ~ {pairs$exposure[i]} + {pairs$Components[i]} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri")

  # fit model for main dataset
  fit2 = glm(as.formula(coding2), data = data, family = binomial(link = "logit"))

  # summary(fit1)
  # summary(fit2)

  medi = mediate(fit1, fit2, treat = glue("{pairs$exposure[i]}"), mediator = glue("{pairs$Components[i]}"), sims = 5000)
  summary(medi)
  
  IQR = quantile(as.matrix(data[[glue("{pairs$exposure[i]}")]]),0.75) - quantile(as.matrix(data[[glue("{pairs$exposure[i]}")]]),0.25)

  tmp <- tibble::as_tibble(dplyr::bind_cols(
      Exposure = pairs$exposure[i],
      Mediator = pairs$Components[i],
      Mode     = pairs$Mode[i],
      Outcome  = pairs$disease[i],
      mz       = pairs$mz[i],
      rt       = pairs$time[i],
      Metabolite_name = pairs$Metabolite[i],
      medi$tau.coef*IQR, medi$tau.ci[1]*IQR,   medi$tau.ci[2]*IQR,   medi$tau.p,  # Total effect
      medi$d.avg*IQR,    medi$d.avg.ci[1]*IQR, medi$d.avg.ci[2]*IQR, medi$d.avg.p,# Indirect effect
      medi$z.avg*IQR,    medi$z.avg.ci[1]*IQR, medi$z.avg.ci[2]*IQR, medi$z.avg.p,# Direct effect
      medi$n.avg*IQR,    medi$n.avg.ci[1]*IQR, medi$n.avg.ci[2]*IQR, medi$n.avg.p # PropMediation
    ))
  colnames(tmp) <- c("Exposure", "Mediator", "Mode", "Outcome", "mz", "rt", "Metabolite_name",
                     "TE",  "TE.LCI",  "TE.UCI",  "TE.Pvalue",
                     "IDE", "IDE.LCI", "IDE.UCI", "IDE.Pvalue",
                     "DE",  "DE.LCI",  "DE.UCI",  "DE.Pvalue",
                     "PropMedi", "PropMedi.LCI", "PropMedi.UCI", "PropMedi.Pvalue")
  return(tmp)
}

# save workspace
save.image(str_c(outpath, "/!Data_for_mediation.rdata"))

load(str_c(outpath, "/!Data_for_mediation.rdata"))



# execute 
tictoc::tic()
safe_pw <- purrr::possibly(Pairwise_1, otherwise = NULL, quiet = FALSE)
res.medi.1 <- purrr::map_dfr(1:nrow(pairs), ~{
  out <- safe_pw(.x)
  if (is.null(out)) message("Pair ", .x, " failed - skipped")
  out
})
tictoc::toc() # ~40 min


# export
res.medi.1 %>%
  fwrite(str_c(outpath, "/DataAnalysis_8_Mediation_analysis",
              "/1_mediation_results.csv"))

# ###########below are code for EMERGE server ####################
# systemd-run --user \
#   --unit="mediation" \
#   --working-directory="/home/mren30/Brain_AD" \
#   -p AllowedCPUs=0-15 \
#   -p MemoryHigh=infinity \
#   -p MemoryMax=infinity \
#   bash -c "bash medi.sh > log.1_mediation_re.log 2>&1"
# systemctl --user status mediation


medi <- fread(str_c(outpath, "/DataAnalysis_8_Mediation_analysis",
              "/1_mediation_results_revised.csv"))

medi %>%
  group_by(Mediator, Mode, Outcome) %>% 
  dplyr::mutate(IDE_fdr = p.adjust(IDE.Pvalue, method = "fdr")) %>%
  ungroup() %>% 
  fwrite(str_c(outpath, "/DataAnalysis_8_Mediation_analysis",
              "/1_mediation_results_revised_new.csv"))

