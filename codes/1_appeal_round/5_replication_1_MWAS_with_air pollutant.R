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
fread(str_c(outpath, "/Met2025_C18_for.MWAS.csv")) -> c18
c18[, -1] %>% 
  dplyr::mutate(across(all_of(matches("V\\d+")), ~ifelse(is.na(.x), min(.x/2, na.rm = TRUE), .x))) -> c18

fread(str_c(outpath, "/Met2025_HILIC_for.MWAS.csv")) -> hilic
hilic[, -1] %>% 
  dplyr::mutate(across(all_of(matches("V\\d+")), ~ifelse(is.na(.x), min(.x/2, na.rm = TRUE), .x))) -> hilic

# merge data
demo %>%
  dplyr::mutate(GUID = as.character(GUID)) %>%
  left_join(ap.data %>% dplyr::select(-GUID), by = "Case_Number") %>% 
  left_join(c18, by = "GUID") %>% 
  as_tibble() -> df.c18

demo %>% 
  dplyr::mutate(GUID = as.character(GUID)) %>%
  left_join(ap.data %>% dplyr::select(-GUID), by = "Case_Number") %>% 
  left_join(hilic, by = "GUID") %>% 
  as_tibble() -> df.hilic


### remove those have been used in main analysis
# import those have been analyzed in discovery phase
df.dis = fread(str_c(work.dir, "/3_Outputs/Cleaned_TRAP_data_for_brain_donor_N162.csv"))

# GUID in discovery cohort
df.dis %>% 
  dplyr::filter(!is.na(SEDAC_PM25NO3_1y)) %>% 
  .$GUID -> guid.dis

df.c18 %>%
  dplyr::filter(!is.na(V1)) %>% # those have metabolomics data -- 347 in total
  dplyr::filter(!is.na(SEDAC_PM25NO3_1y)) %>% # those have air pollution data -- 300 in total
  dplyr::filter(!(GUID %in% guid.dis)) -> df.c18.vali # those also was not used in main analysis -- 143 in total

df.hilic %>%
  dplyr::filter(!is.na(V1)) %>% # those have metabolomics data -- 347 in total
  dplyr::filter(!is.na(SEDAC_PM25NO3_1y)) %>% # those have air pollution data -- 300 in total
  dplyr::filter(!(GUID %in% guid.dis)) -> df.hilic.vali # those also was not used in main analysis -- 143 in total
### remove those have been used in main analysis end 

# save workspace
save.image(str_c(outpath, "/!Data_for_data_analysis_1.rdata"))

###You are now ready for the MWAS data analysis----
load(str_c(outpath, "/!Data_for_data_analysis_1.rdata"))


# 2- Metabolomics-wide Association Study (MWAS)------------------------------
Ex.coef = function(
  expo = "SEDAC_PM25_1y" # air pollution features
  ){
  
  data_c18 = df.c18.vali %>% 
    dplyr::filter(!is.na(df.c18.vali[[expo]]) & !is.na(V1))
  data_hilic = df.hilic.vali %>% 
    dplyr::filter(!is.na(df.hilic.vali[[expo]]) & !is.na(V1))
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


# 2.1.2 for all SEDAC exposures -----------------------------------------
df.c18.vali %>% 
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
saveRDS(out.TRAPtoM.SEDAC, file = str_c(outpath, "/AP_MWAS", 
                                   "/MWAS_TRAP_Metabolome_SEDAC.rds"))

# Save image
save.image(str_c(outpath, "/AP_MWAS", 
                          "/!Data_for_DataAnalysis_3.rdata"))

#------------------------------------------------------------------
# mixture MWAS
load(str_c(outpath, "/AP_MWAS", "/!Data_for_DataAnalysis_3.rdata"))

# SEDAC exposures in mixture, take c18 as example
df.c18 %>% 
  dplyr::select(contains("SEDAC") & !matches("_15y$")) %>% colnames() -> sedac.list
sedac.list %>% 
  subset(str_detect(., "_NO2_|_O3_|_PM25_")) -> sedac.main.components.list
sedac.list %>% 
  subset(str_detect(., "_PM25NH4_|_PM25NO3_|_PM25OC_|_PM25EC_|_PM25SO4_")) -> sedac.pm25.components.list
sedac.list %>% 
  subset(str_detect(., "_PM2.5_")) -> sedac.metal.components.list 
c(sedac.pm25.components.list, sedac.metal.components.list) -> sedac.pm25.full.components.list

# create R6 object for storing results
res <- R6::R6Class(
  "res",
  lock_class = FALSE,
  lock_objects = FALSE
)

# function for mixture analysis: first visit-------------------------
# test = df.c18[, 1:152]
# colnames(test)

Ex.qg.noboot = function(window = "10y", 
                        mix.group = "sedac.main"){                     
  # assign data
  data_c18 = df.c18.vali %>% dplyr::filter(!is.na(V1))  # df.c18
  data_hilic = df.hilic.vali %>% dplyr::filter(!is.na(V1)) # df.hilic

  #loop function
  # mname = "V2"
  F_1 = function(mname){
    # assign equation
    switch(mix.group, 
      "sedac.main"={
        coding = glue("{mname} ~ SEDAC_NO2_{window} + SEDAC_O3_{window} + SEDAC_PM25_{window} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri + death_yr")
        mix.list = sedac.main.components.list  %>% subset(str_detect(., window))
      }, 
      "sedac.pm25"={
        coding = glue("{mname} ~ SEDAC_PM25NH4_{window} + SEDAC_PM25NO3_{window} + SEDAC_PM25OC_{window} + SEDAC_PM25EC_{window} + SEDAC_PM25SO4_{window} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri + death_yr")
        mix.list = sedac.pm25.components.list  %>% subset(str_detect(., window))
      },
      "sedac.metal"={
        coding = glue("{mname} ~ SEDAC_PM2.5_Bromine_{window} + SEDAC_PM2.5_Calcium_{window} + SEDAC_PM2.5_Copper_{window} + SEDAC_PM2.5_Iron_{window} + SEDAC_PM2.5_Lead_{window} + SEDAC_PM2.5_Nickel_{window} + SEDAC_PM2.5_Potassium_{window} + SEDAC_PM2.5_Silicon_{window} + SEDAC_PM2.5_Vanadium_{window} + SEDAC_PM2.5_Zinc_{window} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri + death_yr")
        mix.list = sedac.metal.components.list  %>% subset(str_detect(., window))
      }, 
      "sedac.pm25.full"={
        coding = glue("{mname} ~ SEDAC_PM25NH4_{window} + SEDAC_PM25NO3_{window} + SEDAC_PM25OC_{window} + SEDAC_PM25EC_{window} + SEDAC_PM25SO4_{window} +SEDAC_PM2.5_Bromine_{window} + SEDAC_PM2.5_Calcium_{window} + SEDAC_PM2.5_Copper_{window} + SEDAC_PM2.5_Iron_{window} + SEDAC_PM2.5_Lead_{window} + SEDAC_PM2.5_Nickel_{window} + SEDAC_PM2.5_Potassium_{window} + SEDAC_PM2.5_Silicon_{window} + SEDAC_PM2.5_Vanadium_{window} + SEDAC_PM2.5_Zinc_{window} + Age_MY + sex + race + edu_cat + PMI_MY + ADI_NATRANK + ApoE_tri + death_yr")
        mix.list = sedac.metal.components.list  %>% subset(str_detect(., window))

      })
    
    # qgcomp
    set.seed(1234) # for reproducibility
    qgfit <- qgcomp.glm.noboot(
      as.formula(coding),
      expnms = mix.list,
      data = df,
      family = gaussian(),  # linear model
      q = 4  # Quartiles for exposure quantiles
    )
    # summary(qgfit)
    #extract coefficients
    weights = qgfit$fit$coefficients %>% .[names(.) %in% mix.list]

    c(qgfit$psi, qgfit$ci, sum(weights[which(weights>=0)]), sum(weights[which(weights<=0)]), qgfit$tstat[2], qgfit$pval[2]) %>% 
          t() %>% 
          as_tibble() %>% 
          dplyr::mutate(Windows = window,
                        Outcome = mname) %>% 
          purrr::set_names("MixEffect", "LowerCI", "UpperCI", 
                           "Positive.Size", "Negative.Size", "tstat", "P.value", "Windows", "Outcome") %>% 
          dplyr::select(Outcome, Windows, everything()) -> qg.coef
    # extract weights
    if(qg.coef$Positive.Size[[1]]!=0){
      pos.weights = weights[which(weights>=0)]/sum(weights[which(weights>=0)])
      pos.weights %>%
          as_tibble() %>% cbind(names(pos.weights)) %>% 
          dplyr::mutate(Direction = "Positive") %>% 
          purrr::set_names("Weight", "Component", "Direction") -> qg.pos
    }else{qg.pos=tibble(Component = NA, Weight = NA, Direction = NA)}
    if(qg.coef$Negative.Size[[1]]!=0){
      neg.weights = weights[which(weights<=0)]/sum(weights[which(weights<=0)])
      neg.weights %>%
          as_tibble() %>% cbind(names(neg.weights)) %>% 
          dplyr::mutate(Direction = "Negative") %>% 
          purrr::set_names("Weight", "Component", "Direction") -> qg.neg
    }else{qg.neg=tibble(Component = NA, Weight = NA, Direction = NA)}
    # rbind
    rbind(qg.pos, qg.neg) %>% 
        dplyr::mutate(Windows = window,
                      Outcome = mname) %>% 
        dplyr::select(Outcome, Windows, Component, Direction, Weight) -> qg.weight
    return(list(qg.coef, qg.weight))
  }# end F_1
  
  # c18 ---------------
  data_c18 %>% dplyr::select(matches("V\\d+")) -> M
  data_c18 %>% dplyr::select(contains(glue("_{window}")))  -> X
  data_c18 %>% dplyr::select(all_of(c("Age_MY",  "sex", "race", "edu_cat", "PMI_MY", "ADI_NATRANK", "ApoE_tri", "death_yr"))) -> C
  df <- dplyr::bind_cols(X, C, M)
  mname <- colnames(M) 
  # execute for c18
  furrr::future_map(mname, ~F_1(.x)) -> tmp
  print(glue("MWAS for the relations between [ {mix.group} ] Air pollution mixture of [ {window} ] window and metabolome(c18) have been done!"))
  # store coefficients
  purrr::map(tmp, 1) %>% 
    bind_rows() -> res[[mix.group]]$c18$coef[[window]]
  # store weights
  purrr::map(tmp, 2) %>% 
    bind_rows() -> res[[mix.group]]$c18$weight[[window]]    

  # hilic ---------------
  data_hilic %>% dplyr::select(matches("V\\d+")) -> M
  data_hilic %>% dplyr::select(contains(glue("_{window}")))  -> X
  data_hilic %>% dplyr::select(all_of(c("Age_MY",  "sex", "race", "edu_cat", "PMI_MY", "ADI_NATRANK", "ApoE_tri", "death_yr"))) -> C
  df <- dplyr::bind_cols(X, C, M)
  mname <- colnames(M) 
  # execute for hilic
  furrr::future_map(mname, ~F_1(.x)) -> tmp
  print(glue("MWAS for the relations between [ {mix.group} ] Air pollution mixture of [ {window} ] window and metabolome(hilic) have been done!"))
  # store coefficients
  purrr::map(tmp, 1) %>% 
    bind_rows() -> res[[mix.group]]$hilic$coef[[window]]
  # store weights
  purrr::map(tmp, 2) %>% 
    bind_rows() -> res[[mix.group]]$hilic$weight[[window]]    
  return(res)
}

# execute--------
windows = c("1y", "3y", "5y", "10y")
mix.groups = c("sedac.main", "sedac.pm25", "sedac.metal", "sedac.pm25.full")

# execute
tictoc::tic()
purrr::map2(rep(windows, each = length(mix.groups)),
            rep(mix.groups, time = length(windows)),
            Ex.qg.noboot)
tictoc::toc() # xxx

# save rds file for further extraction
saveRDS(res, file = str_c(outpath, "/AP_MWAS", 
                                   "/MWAS_TRAP_Metabolome_mixture.rds"))

# Save image
save.image(str_c(outpath, "/AP_MWAS", 
                          "/!Data_for_DataAnalysis_4.rdata"))






##############################
# import MWAS results
out.TRAPtoM.SEDAC = readRDS(str_c(outpath, "/AP_MWAS", 
                          "/MWAS_TRAP_Metabolome_SEDAC.rds"))
out.TRAPtoM.mixture = readRDS(str_c(outpath, "/AP_MWAS", 
                          "/MWAS_TRAP_Metabolome_mixture.rds"))
# SEDAC----
# c18
map2(out.TRAPtoM.SEDAC, out.TRAPtoM.SEDAC %>% list.names(), 
     ~.x$c18 %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = .y) %>% 
      dplyr::select(Exposure, P_fdr, Pvalue_a, everything())) %>% 
rbindlist()  -> c18.SEDAC

# hilic
map2(out.TRAPtoM.SEDAC, out.TRAPtoM.SEDAC %>% list.names(), 
     ~.x$hilic %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = .y) %>% 
      dplyr::select(Exposure, P_fdr, Pvalue_a, everything())) %>% 
rbindlist()  -> hilic.SEDAC

# extract p value counts
map_dfr(out.TRAPtoM.SEDAC %>% list.names(), 
  ~tibble(
    Exposure = .x,  # should be performed by various phenotypes here!!
    Group = "c18", 
    Type = "SEDAC",
    P_raw_0.05 =  sum(c18.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.05, na.rm = T),
    P_raw_0.01 =  sum(c18.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.01, na.rm = T),
    P_raw_0.001 = sum(c18.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(c18.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(c18.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(c18.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ) %>% 
  rbind(
  tibble(
    Exposure= .x, 
    Group = "hilic", 
    Type = "SEDAC",
    P_raw_0.05 =  sum(hilic.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.05, na.rm = T),
    P_raw_0.01 =  sum(hilic.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.01, na.rm = T),
    P_raw_0.001 = sum(hilic.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$Pvalue_a<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(hilic.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(hilic.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(hilic.SEDAC %>% dplyr::filter(Exposure == .x) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ))) %>% 
fwrite(str_c(outpath, "/AP_MWAS", 
               "/Pvalue_counts_MWAS.csv"))



#################################
# mixtures
mix.groups = c("sedac.main", "sedac.pm25", "sedac.metal", "sedac.pm25.full")
windows = c("1y", "3y", "5y", "10y")
## SEDAC main----
# c18
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.main"]][["c18"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())
            }
          ) %>% 
            rbindlist()  -> c18.sedac.main
# hilic
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.main"]][["hilic"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())
            }
          ) %>% 
            rbindlist()  -> hilic.sedac.main
# extract p value counts
map_dfr(windows, 
  ~{tibble(
    Exposure = str_c("mix_", .x), 
    Group = "c18", 
    Type = "SEDAC.main",
    P_raw_0.05 =  sum(c18.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(c18.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(c18.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(c18.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(c18.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(c18.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ) %>% 
  rbind(
  tibble(
    Exposure= str_c("mix_", .x), 
    Group = "hilic", 
    Type = "SEDAC.main",
    P_raw_0.05 =  sum(hilic.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(hilic.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(hilic.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(hilic.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(hilic.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(hilic.sedac.main %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ))
  }
) %>% 
fwrite(str_c(outpath, "/AP_MWAS", 
               "/Pvalue_counts_MWAS_mixture_main pollutant.csv"))

## SEDAC PM2.5 components----
# c18
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.pm25"]][["c18"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())}
          ) %>% 
            rbindlist()  -> c18.sedac.pm25
# hilic
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.pm25"]][["hilic"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())
            }
          ) %>% 
            rbindlist()  -> hilic.sedac.pm25
# extract p value counts
map_dfr(windows, 
  ~{tibble(
    Exposure = str_c("mix_", .x), 
    Group = "c18", 
    Type = "SEDAC.pm25",
    P_raw_0.05 =  sum(c18.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(c18.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(c18.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(c18.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(c18.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(c18.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ) %>% 
  rbind(
  tibble(
    Exposure= str_c("mix_", .x), 
    Group = "hilic", 
    Type = "SEDAC.pm25",
    P_raw_0.05 =  sum(hilic.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(hilic.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(hilic.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(hilic.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(hilic.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(hilic.sedac.pm25 %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ))
  }
) %>% 
fwrite(str_c(outpath, "/AP_MWAS", 
               "/Pvalue_counts_MWAS_mixture_pm25 species.csv"))

## SEDAC metal components----
# c18
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.metal"]][["c18"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())
            }
          ) %>% 
            rbindlist()  -> c18.sedac.metal
# hilic
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.metal"]][["hilic"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())
            }
          ) %>% 
            rbindlist()  -> hilic.sedac.metal
# extract p value counts
map_dfr(windows, 
  ~{tibble(
    Exposure = str_c("mix_", .x), 
    Group = "c18", 
    Type = "SEDAC.metal",
    P_raw_0.05 =  sum(c18.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(c18.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(c18.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(c18.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(c18.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(c18.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ) %>% 
  rbind(
  tibble(
    Exposure= str_c("mix_", .x), 
    Group = "hilic", 
    Type = "SEDAC.metal",
    P_raw_0.05 =  sum(hilic.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(hilic.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(hilic.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(hilic.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(hilic.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(hilic.sedac.metal %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted p
  ))
  }
) %>% 
fwrite(str_c(outpath, "/AP_MWAS", 
               "/Pvalue_counts_MWAS_mixture_pm25 metal.csv"))

# pm2.5 full components
# c18
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.pm25.full"]][["c18"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())
            }
          ) %>% 
            rbindlist()  -> c18.sedac.pm25.full
# hilic
purrr::map(windows, 
          ~ {out.TRAPtoM.mixture[["sedac.pm25.full"]][["hilic"]][["coef"]][[.x]] %>% 
              dplyr::mutate(P_fdr = p.adjust(P.value, method="BH"), 
                            Exposure = str_c("mix_", .x)) %>% 
              dplyr::select(Exposure, P_fdr, P.value, everything())
            }
          ) %>% 
            rbindlist()  -> hilic.sedac.pm25.full
# extract p value counts
map_dfr(windows, 
  ~{tibble(
    Exposure = str_c("mix_", .x), 
    Group = "c18", 
    Type = "SEDAC.pm25.full",
    P_raw_0.05 =  sum(c18.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(c18.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(c18.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(c18.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(c18.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(c18.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted q
  ) %>% 
  rbind(
  tibble(
    Exposure= str_c("mix_", .x), 
    Group = "hilic", 
    Type = "SEDAC.pm25.full",
    P_raw_0.05 =  sum(hilic.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.05, na.rm = T),
    P_raw_0.01 =  sum(hilic.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.01, na.rm = T),
    P_raw_0.001 = sum(hilic.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P.value<0.001, na.rm = T), # unadjusted p
    P_fdr_0.2 =   sum(hilic.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.2, na.rm = T),
    P_fdr_0.05 =  sum(hilic.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.05, na.rm = T),
    P_fdr_0.01 =  sum(hilic.sedac.pm25.full %>% dplyr::filter(Exposure == str_c("mix_", .x)) %>% .$P_fdr<0.01, na.rm = T) # fdr-adjusted q
  ))
  }
) %>% 
fwrite(str_c(outpath, "/AP_MWAS", 
               "/Pvalue_counts_MWAS_mixture_pm25 full components.csv"))

