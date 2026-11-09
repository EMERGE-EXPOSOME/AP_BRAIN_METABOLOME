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
outpath <- str_c(work.dir, "/3_Outputs")


# Create output folder for the upcoming MWAS results of traffic-related air pollution and metabolomic signatures
if(!file.exists(str_c(outpath,"/DataAnalysis_1_TRAP_Metabolome"))){
  dir.create(str_c(outpath,"/DataAnalysis_1_TRAP_Metabolome"))
}

# 2-Load data ----------
# load demographic data
fread(str_c(outpath, "/Cleaned_data_demographics.csv")) -> demo
demo %<>% 
  dplyr::mutate(ApoE_tri = if_else(ApoE_tri == "", NA, ApoE_tri))
demo$Case_Number
demo$date_of_death

# load air pollution data
fread(str_c(outpath, "/Cleaned_data_Air_pollution(CACES).csv")) -> ap.data

# load metabolomic data
fread(str_c(outpath, "/Cleaned_data_for_MWAS_c18.csv")) %>% 
  dplyr::mutate(Case_Number = str_replace(Sample_ID, "_", "-")) -> c18
c18 %<>% 
  dplyr::mutate(across(all_of(matches("V\\d+")), ~ifelse(is.na(.x), min(.x/2, na.rm = TRUE), .x)))
fread(str_c(outpath, "/Cleaned_data_for_MWAS_hilic.csv")) %>% 
  dplyr::mutate(Case_Number = str_replace(Sample_ID, "_", "-")) -> hilic
hilic %<>% 
  dplyr::mutate(across(all_of(matches("V\\d+")), ~ifelse(is.na(.x), min(.x/2, na.rm = TRUE), .x)))

c18$Case_Number %>% str_detect("^q3") %>% sum()

# merge data
demo %>%
  left_join(ap.data[, -9], by = "Case_Number") %>% 
  left_join(c18, by = "Case_Number") %>% 
  as_tibble() -> df.c18

demo %>% 
  left_join(ap.data[, -9], by = "Case_Number") %>% 
  left_join(hilic, by = "Case_Number") %>% 
  as_tibble() -> df.hilic

save.image(str_c(outpath, "/!Data_for_data_analysis_2(CACES).rdata"))

###You are now ready for the MWAS data analysis----
load(str_c(outpath, "/!Data_for_data_analysis_2(CACES).rdata"))

# df.c18 %>% 
#   dplyr::filter(!is.na(df.c18[["Rline_PM25_1y"]]) & !is.na(V1)) %>% nrow()
df.c18 %>% nrow()
colnames(df.c18)[1:150]
df.c18 %>% 
  dplyr::filter(!is.na(V1)) %>% 
  dplyr::select(Case_Number, OBJECTID.x, GUID.x, CRIN_number, date_of_death.x:pm10_10y) %>% 
  dplyr::rename(GUID = GUID.x, 
                OBJECTID = OBJECTID.x) -> out
out$dist_km %>% hist()

out %>% 
  fwrite(str_c(outpath, "/Cleaned_TRAP_data_for_brain_donor_N162(CACES).csv"))


# 2- Metabolomics-wide Association Study (MWAS)------------------------------

Ex.coef = function(
  expo = "pm10_1y" # air pollution features
  ){
  
  data_c18 = df.c18 %>% 
    dplyr::filter(!is.na(df.c18[[expo]]) & !is.na(V1))
  data_hilic = df.hilic %>% 
    dplyr::filter(!is.na(df.c18[[expo]]) & !is.na(V1))
  print(glue("For TRAP [ {expo} ], a total of [ {nrow(data_c18)} ] obs were obtained for downstream analysis!"))
  
  #create list to store results
  res = list(l1=list(), l2 = list())
  names(res) = c("c18", "hilic")

  # c18------------
  data_c18 %>% dplyr::select(all_of(expo)) -> X
  data_c18 %>% dplyr::select(matches("V\\d+")) -> M
  data_c18 %>% dplyr::select(Case_Number:death_mt) -> C
  df <- dplyr::bind_cols(X, M, C)
  name <- colnames(M)

  #loop function
  F1 = function(
      meta,# metabolites
      name # name of metabolites
  ){
    tmp = NULL
    # name = "V6465"
    # meta = df$V6465

      #continuous exposures, IQR change
      IQR = quantile(as.matrix(df[,expo]),0.75) - quantile(as.matrix(df[,expo]),0.25)
      #1. crude Coef (cCoef) and 95% CI
      eqt = paste(name, " ~ ", expo)
      fit <- lm(as.formula(eqt), data = df)
      tmp <- summary(fit)$coefficients
      if(is.null(tmp) | !(expo %in% row.names(tmp))){
        # in some cases, model may not fit (due to low variance of specific metabolites)
        # so we need to assign tmp1 as follow:
        tmp1 <- cbind(cCoef = NA,
                      LCI_c = NA,
                      UCI_c = NA,
                      t_c = NA,
                      Pvalue_c = NA)
        print(str_c("Warning!! metabolite [", name, "] could not fit raw models for ", expo, " due to low variates !"))
      }else{
        tmp1 <- cbind(cCoef = coef(fit)[[2]]*IQR,
                      LCI_c = confint.default(fit)[2,1]*IQR,
                      UCI_c = confint.default(fit)[2,2]*IQR, 
                      t_c = tmp[2,3],
                      Pvalue_c = tmp[2,4])
      }

      #2. adjusted aCoef and 95% CI
      eqt = paste(name, " ~ ", expo, " + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri + death_yr")

      fit2 <- lm(as.formula(eqt), data = df)
      tmp <- summary(fit2)$coefficients
      if(is.null(tmp) | !(expo %in% row.names(tmp))){
        # in some cases, model may not fit (due to low variance of specific metabolites)
        # so we need to assign tmp1 as follow:
        tmp2 <- cbind(aCoef = NA,
                      LCI_a = NA,
                      UCI_a = NA, 
                      t_a = NA,
                      Pvalue_a = NA)
        print(str_c("Warning!! metabolite [", name, "] could not fit adjusted models for ", expo, " due to low variates !"))
      }else{      
        tmp2 <- cbind(aCoef = coef(fit2)[[2]]*IQR,
                      LCI_a = confint.default(fit2)[2,1]*IQR,
                      UCI_a = confint.default(fit2)[2,2]*IQR, 
                      t_a = tmp[2,3],
                      Pvalue_a = tmp[2,4])
      }
      #3. combine cCoef and aCoef
      cbind(tmp1, tmp2) %>% 
          tibble::as_tibble() %>% 
          dplyr::mutate(Components = name) %>% 
          dplyr::select(Components, everything()) -> tmp3
      return(tmp3)
    }

  # A TEST
  # purrr::map2_dfr(M[,6465:6470],colnames(M[,6465:6470]), ~F1(.x, .y))
  
  # execute for c18 column 
  res[["c18"]] <- purrr::map2_dfr(M, name, F1)
  print(paste("Associatios between [c18] metabolites and ", expo," on  have been calculated!"))
  gc()

  # hilic----------------
  data_hilic %>% dplyr::select(all_of(expo)) -> X
  data_hilic %>% dplyr::select(matches("V\\d+")) -> M
  data_hilic %>% dplyr::select(Case_Number:death_mt) -> C
  df <- dplyr::bind_cols(X, M, C)
  name <- colnames(M)

  # A TEST
  # purrr::map2_dfr(M[,6465:6470],colnames(M[,6465:6470]), ~F1(.x, .y))
  # execute for hilic column 
  res[["hilic"]] <- purrr::map2_dfr(M, name, F1)
  print(paste("Associatios between [hilic] metabolites and ", expo," on  have been calculated!"))
  gc()

  return(res)
}

# 2.1.1 for all CACES exposures ------------------------------------------
df.c18 %>% 
  dplyr::select(contains("pm25") | contains("pm10")) %>% colnames() -> caces.list

# execute
tictoc::tic()
plan(multisession, workers = 8)
furrr::future_map(.x = caces.list, 
                  ~Ex.coef(.x)) -> out.TRAPtoM
names(out.TRAPtoM) = caces.list
plan(sequential)
tictoc::toc() # ~ 1 hours
gc()

# save rds file for further extraction
saveRDS(out.TRAPtoM, file = str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                                   "/MWAS_TRAP_Metabolome_CACES.rds"))

# Save image
save.image(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                          "/!Data_for_DataAnalysis_2(CACES).rdata"))




# 3- Export counts of pvalues for each set of XtoM----------------
# load rds file
res <- readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                                   "/MWAS_TRAP_Metabolome_CACES.rds"))

# c18
map2(res, res %>% list.names(), 
     ~.x$c18 %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = .y) %>% 
      dplyr::select(Exposure, P_fdr, Pvalue_a, everything())) %>% 
rbindlist()  -> c18.CACES
# hilic
map2(res, res %>% list.names(), 
     ~.x$hilic %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = .y) %>% 
      dplyr::select(Exposure, P_fdr, Pvalue_a, everything())) %>% 
rbindlist()  -> hilic.CACES

# extract p value counts
map_dfr(res %>% list.names(), 
  ~tibble(
    Exposure = .x, 
    Group = "c18", 
    Type = "CACES",
    P_raw_0.05 =  sum(c18.CACES %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.05, na.rm = T),
    P_raw_0.01 =  sum(c18.CACES %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.01, na.rm = T),
    P_raw_0.001 = sum(c18.CACES %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(c18.CACES %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(c18.CACES %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(c18.CACES %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ) %>% 
  rbind(
  tibble(
    Exposure= .x, 
    Group = "hilic", 
    Type = "CACES",
    P_raw_0.05 =  sum(hilic.CACES %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.05, na.rm = T),
    P_raw_0.01 =  sum(hilic.CACES %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.01, na.rm = T),
    P_raw_0.001 = sum(hilic.CACES %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(hilic.CACES %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(hilic.CACES %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(hilic.CACES %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ))) %>% 
fwrite(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
               "/Pvalue_counts_XtoM_CACES.csv"))
