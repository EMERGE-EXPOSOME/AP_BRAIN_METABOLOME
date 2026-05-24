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
# outpath <- str_c(getwd(), "")

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

# load air pollution data
fread(str_c(outpath, "/Cleaned_data_Air_pollution.csv")) -> ap.data

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
  left_join(ap.data, by = "Case_Number") %>% 
  left_join(c18, by = "Case_Number") %>% 
  as_tibble() -> df.c18

demo %>% 
  left_join(ap.data, by = "Case_Number") %>% 
  left_join(hilic, by = "Case_Number") %>% 
  as_tibble() -> df.hilic

save.image(str_c(outpath, "/!Data_for_data_analysis_1.rdata"))
###You are now ready for the MWAS data analysis----
load(str_c(outpath, "/!Data_for_data_analysis_1.rdata"))

# df.c18 %>% 
#   dplyr::filter(!is.na(df.c18[["Rline_PM25_1y"]]) & !is.na(V1)) %>% nrow()
df.c18 %>% nrow()
colnames(df.c18)[1:150]
df.c18 %>% 
  dplyr::filter( !is.na(V1))  %>% 
  dplyr::select(Case_Number,OBJECTID.x, GUID.x, CRIN_number, date_of_death:SEDAC_PM2.5_Zinc_15y) %>% 
  dplyr::rename(GUID = GUID.x, 
                OBJECTID = OBJECTID.x) -> out
out$death_yr %>% range(na.rm = T)
out %>% 
  fwrite(str_c(outpath, "/Cleaned_TRAP_data_for_brain_donor_N161.csv"))


df.c18 %>% 
  dplyr::filter(!is.na(df.c18[["SEDAC_PM25NO3_1y"]]) & !is.na(V1)) -> tmp

table(tmp$ApoE_tri, tmp$Braak_Stage_new)
table(tmp$ApoE_tri, tmp$ABC_new)
table(tmp$ApoE_tri, tmp$CERAD_new)
# of note: in 2.3.1, here we did not exclude participants without TRAP data
tt = c("Age_MY", "sex", "race", "edu_cat", "PMI_MY",
       "ADI_NATRANK", "ApoE_tri", "ABC_new", 
       "Braak_Stage_new", "CERAD_new") 
cat = c("sex", "race", "edu_cat", 
       "ApoE_tri", "Braak_Stage_new") 

# create TableOne: AllPTBcsAllTerm
tab1 <-  print(CreateTableOne(vars = tt, # here we use All term birth as control group, to show whole participants in our study
                              data = tmp, factorVars = cat,
                              includeNA = FALSE),
               showAllLevels = TRUE)
rownames(tab1)  %>% as_tibble()  %>%
  cbind(tab1 %>% as_tibble()) %>%
  as_tibble() %>%
  set_names(c("Variables", "Levels","Overall")) %>% 
writexl::write_xlsx(stringr::str_c(outpath,
                                  "/TableOne_1_Descriptive_results_N161.xlsx"))


# 2- Metabolomics-wide Association Study (MWAS)------------------------------

Ex.coef = function(
  expo = "SEDAC_PM25_5y" # air pollution features
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

# 2.1.1 for all SEDAC exposures -----------------------------------------
df.c18 %>% 
  dplyr::select(contains("SEDAC") & !matches("_15y$")) %>% colnames() -> SEDAC.list
# execute
tictoc::tic()
plan(multisession, workers = 8)
furrr::future_map(.x = SEDAC.list, 
                  ~Ex.coef(.x)) -> out.TRAPtoM.SEDAC
names(out.TRAPtoM.SEDAC) = SEDAC.list
plan(sequential)
tictoc::toc() # ~3 hours
gc()
# save rds file for further extraction
saveRDS(out.TRAPtoM.SEDAC, file = str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                                   "/MWAS_TRAP_Metabolome_SEDAC.rds"))


# Save image
save.image(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                          "/!Data_for_DataAnalysis_2.rdata"))


