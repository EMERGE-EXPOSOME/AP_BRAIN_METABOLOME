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
load(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                          "/!Data_for_DataAnalysis_2.rdata"))
ls()
###You are now ready for the MWAS data analysis----
colnames(df.c18)[1:150]

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

Ex.qg.noboot = function(window = "3y", 
                        mix.group = "sedac"){                     
  # assign data
  data_c18 = df.c18 %>% dplyr::filter(!is.na(V1))  # df.c18
  data_hilic = df.hilic %>% dplyr::filter(!is.na(V1)) # df.hilic

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
    }else{qg.pos=NULL}
    if(qg.coef$Negative.Size[[1]]!=0){
      neg.weights = weights[which(weights<=0)]/sum(weights[which(weights<=0)])
      neg.weights %>%
          as_tibble() %>% cbind(names(neg.weights)) %>% 
          dplyr::mutate(Direction = "Negative") %>% 
          purrr::set_names("Weight", "Component", "Direction") -> qg.neg
    }else{qg.neg=NULL}
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
  map(tmp, 1) %>% 
    bind_rows() -> res[[mix.group]]$c18$coef[[window]]
  # store weights
  map(tmp, 2) %>% 
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
  map(tmp, 1) %>% 
    bind_rows() -> res[[mix.group]]$hilic$coef[[window]]
  # store weights
  map(tmp, 2) %>% 
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
saveRDS(res, file = str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                                   "/MWAS_TRAP_Metabolome_mixture.rds"))



# Save image
save.image(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome", 
                          "/!Data_for_DataAnalysis_2.rdata"))
