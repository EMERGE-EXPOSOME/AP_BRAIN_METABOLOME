rm(list=(ls()))
# 1- Prepare steps -------------------------------------
# load package ---------
pacman::p_load("readxl",
               "writexl",
               "purrr",
               "furrr",
               "lubridate", 
               "dplyr", 
               "magrittr",
               "stringr", 
               "rlist", 
               "tidyverse", 
               "tableone", 
               "gghalves", 
               "ddpcr",
               "glue", 
               "qgcomp", 
               "data.table", 
               "lme4", 
               "lmerTest", 
               "statmod", 
               "R6", 
			         "tictoc", 
               "metapone", 
               "kableExtra", 
               "RColorBrewer", 
               "ggtext")

# set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")
# outpath <- str_c(getwd(), "")

# Create output folder for the upcoming MWAS results of traffic-related air pollution and metabolomic signatures
if(!file.exists(str_c(outpath,"/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/TRAP_mixture"))){
  dir.create(str_c(outpath,"/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/TRAP_mixture"))
}

# load data to be performed using metapone----
out.TRAPtoM.mixture = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
                                  "/MWAS_TRAP_Metabolome_mixture.rds"))
out.TRAPtoM.mixture %>% ls()
out.TRAPtoM.mixture[["sedac.main"]][["c18"]][["coef"]] %>% list.names()
out.TRAPtoM.mixture[["sedac.pm25"]][["c18"]][["coef"]] %>% list.names()
out.TRAPtoM.mixture[["sedac.metal"]][["c18"]][["coef"]] %>% list.names()
out.TRAPtoM.mixture[["sedac.pm25.full"]][["c18"]][["coef"]] %>% list.names()

# 2-metapone for individual air pollution and metabolome----
# prepare pa and hmdbCompMZ datasets for metapone----
pa<-readRDS(str_c(inpath, 
                  "/Metabolomic_data", 
                  "/pa_human.rds"))
# View(pa)
data("hmdbCompMZ")
pa<-pa %>% filter(flag==1) # level-1 list

# set parameters----
pos_adduct <- c("M+2ACN+H", "M+ACN+H", "M+ACN+Na", "M+H", "M+K", "M+Na")
neg_adduct <- c("M-2H", "M-H2O-H", "M-H", "M+Na-2H", "M+Cl")
p_thres <- 0.05  # raw p value into metapone, and we will use adjusted p values in metapone results for display?
n_permu <- 1000

# Initial res to store metapone results----
res_mtp <- R6::R6Class(
  "res",
  lock_class = FALSE,
  lock_objects = FALSE
)


# import untargeted metabolomic mz and retention time information----
fread(str_c(outpath, "/log2.transformed.metabolomics_C18.csv"))  %>% 
    .[, 1:10] -> c18_info
fread(str_c(outpath, "/log2.transformed.metabolomics_HILIC.csv"))  %>% 
    .[, 1:10] -> hilic_info

# FUNCTION to execute metapone for listed files----
F_metapone = function(
  mix.group = "sedac.main", # sedac.main, sedac.pm25, sedac.metal, sedac.pm25.full
  window = "1y"){ # 1y, 3y, 5y, 10y
  
  # import data----
  # c18
  out.TRAPtoM.mixture[[mix.group]][["c18"]][["coef"]][[window]] %>% 
    dplyr::filter(!is.na(Outcome)) %>% 
    dplyr::mutate(mid = as.numeric(str_extract(Outcome, "(?<=V)\\d+")))  -> test
  c18_info %>% 
    dplyr::mutate(mid = row_number()) %>% 
    full_join(test, by = "mid") %>% 
    dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
    dplyr::select(mz, time, P.value, tstat, P_fdr) -> neg
  # hilic
  out.TRAPtoM.mixture[[mix.group]][["hilic"]][["coef"]][[window]] %>% 
    dplyr::filter(!is.na(Outcome)) %>% 
    dplyr::mutate(mid = as.numeric(str_extract(Outcome, "(?<=V)\\d+")))  -> test
  hilic_info %>% 
    dplyr::mutate(mid = row_number()) %>% 
    full_join(test, by = "mid") %>% 
    dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
    dplyr::select(mz, time, P.value, tstat, P_fdr) -> pos

  # data used for metapone
  neg %>% 
    dplyr::select(mz, time, P.value, tstat) %>% 
    set_names(c("mz", "time", "p.value", "t.value")) %>% 
    as.matrix() -> neg
  pos %>% 
    dplyr::select(mz, time, P.value, tstat) %>% 
    set_names(c("mz", "time", "p.value", "t.value")) %>% 
    as.matrix() -> pos
  dat = list(pos, neg)
  mode = list("pos", "neg")
  # metapone execution
  metapone(dat, mode, pa, hmdbCompMZ = hmdbCompMZ,
           pos.adductlist = pos_adduct, neg.adductlist = neg_adduct,
           p.threshold = p_thres, n.permu = n_permu,
           fractional.count.power=0.5, max.match.count=10, use.fgsea = FALSE) -> metapone.fit
  # metapone.fit[[str_c(grp_name, "_v", as.character(visit))]]
  print(str_c("metapone analysis has been done in [", mix.group, "] in [ ", window, " ] !"))
  
  # pathway tables: extract raw-p-value (p_value) <0.05 and significant metabolites (n_significant_metabolites) >1
  ptable(metapone.fit) %>% 
    dplyr::filter(p_value < 0.05 & n_significant_metabolites > 1) %>% 
    arrange(p_value) -> pathways
  # export to xlsx
  pathways %>% 
    dplyr::mutate(Pathways = row.names(pathways)) %>% 
    writexl::write_xlsx(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM", 
                          "/metapone_pathways_mixture_", mix.group, "_", window, ".xlsx"))
  # feature tables: extract features that in each pathways
  ftable(metapone.fit)[row.names(ptable(metapone.fit) %>% dplyr::filter(p_value < 0.05 & n_significant_metabolites > 1))] -> features
  map2_dfr(
    features, 
    features %>% ls(),
    ~ {
      as.data.frame(.x) %>% 
      dplyr::mutate(pathways = as.character(.y)) %>% 
      dplyr::select(pathways, everything())
    }
    ) -> features
  gc()
  return(list(names = str_c(mix.group, "_", window),
              pathways, 
              features))
}

mix.groups = c("sedac.main", "sedac.pm25", "sedac.metal", "sedac.pm25.full")
windows = c("1y", "3y", "5y", "10y")

# execute SEDAC
list2 = mix.groups[2:5]

tictoc::tic()
# plan(multisession, workers = 8)
res.list2 = purrr::map2(rep(list2, each = length(windows)), 
                               rep(windows, times = length(list2)), 
                               F_metapone)
ddpcr::quiet(
purrr::map(res.list2, ~ {
  .x[[2]] -> res_mtp$pathways[[.x[[1]]]]
  .x[[3]] -> res_mtp$features[[.x[[1]]]]
})
)
# plan(sequential) # 10 hours
tictoc::toc()

# save res_mtp
saveRDS(res_mtp, file = str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM", 
               "/Pathway_and_signatures_results_XtoM_mixture.rds"))


# 3- Export counts of pvalues for each set of XtoM----------------
out.TRAPtoM.mixture[["sedac.metal"]][["c18"]][["coef"]][["5y"]] %>% .$P.value %>% hist()
mix.groups
windows

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
fwrite(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/TRAP_mixture", 
               "/Pvalue_counts_XtoM_SEDAC.main_mixture.csv"))

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
fwrite(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/TRAP_mixture", 
               "/Pvalue_counts_XtoM_SEDAC.pm25_mixture.csv"))

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
fwrite(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/TRAP_mixture", 
               "/Pvalue_counts_XtoM_SEDAC.metal_mixture.csv"))

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
fwrite(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_Analyses_using maximum datasets/TRAP_mixture", 
               "/Pvalue_counts_XtoM_SEDAC.pm25.full_mixture.csv"))

