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
if(!file.exists(str_c(outpath,"/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets"))){
  dir.create(str_c(outpath,"/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets"))
}

# load data to be performed using metapone----
out.TRAPtoM.SEDAC = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
                                  "/MWAS_TRAP_Metabolome_SEDAC.rds"))

out.TRAPtoM.SEDAC$SEDAC_PM25_10y
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
  type = "SEDAC",
  grp_name = "CO_1y"){
  
  # assign data
  switch(type, 
    "SEDAC"={
      res = out.TRAPtoM.SEDAC
    })
  # import data----
  # c18
  res[[str_c(type, "_", grp_name)]][["c18"]] %>% 
    dplyr::filter(!is.na(Components)) %>% 
    dplyr::mutate(mid = as.numeric(str_extract(Components, "(?<=V)\\d+")))  -> test
  c18_info %>% 
    dplyr::mutate(mid = row_number()) %>% 
    full_join(test, by = "mid") %>% 
    dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
    dplyr::select(mz, time, Pvalue_a, t_a, P_fdr) -> neg
  # hilic
  res[[str_c(type, "_", grp_name)]][["hilic"]]  %>% 
    dplyr::filter(!is.na(Components)) %>% 
    dplyr::mutate(mid = as.numeric(str_extract(Components, "(?<=V)\\d+"))) -> test
  hilic_info %>% 
    dplyr::mutate(mid = row_number()) %>% 
    full_join(test, by = "mid") %>% 
    dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
    dplyr::select(mz, time, Pvalue_a, t_a, P_fdr) -> pos

  # data used for metapone
  neg %>% 
    dplyr::select(mz, time, Pvalue_a, t_a) %>% 
    set_names(c("mz", "time", "p.value", "t.value")) %>% 
    as.matrix() -> neg
  pos %>% 
    dplyr::select(mz, time, Pvalue_a, t_a) %>% 
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
  print(str_c("metapone analysis has been done in [", grp_name, "] set of ", type, " !"))
  
  # pathway tables: extract raw-p-value (p_value) <0.05 and significant metabolites (n_significant_metabolites) >1
  ptable(metapone.fit) %>% 
    dplyr::filter(p_value < 0.05 & n_significant_metabolites > 1) %>% 
    arrange(p_value) -> pathways
  # export to xlsx
  pathways %>% 
    dplyr::mutate(Pathways = row.names(pathways)) %>% 
    writexl::write_xlsx(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM", 
                          "/metapone_pathways_", type, "_", grp_name, ".xlsx"))
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
  return(list(names = str_c(type, "_", grp_name),
              pathways, 
              features))
}


# execute SEDAC
list2 = out.TRAPtoM.SEDAC %>% list.names() %>% str_extract("(?<=SEDAC_)[a-zA-Z0-9\\_\\.]+")

tictoc::tic()
plan(multisession, workers = 12)
res.list2 = furrr::future_map2(rep("SEDAC", each = length(list2)), 
                              list2, 
                              F_metapone)
ddpcr::quiet(
purrr::map(res.list2, ~ {
  .x[[2]] -> res_mtp$pathways[[.x[[1]]]]
  .x[[3]] -> res_mtp$features[[.x[[1]]]]
})
)
plan(sequential) # xxx hours
tictoc::toc()

# save res_mtp
saveRDS(res_mtp, file = str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM", 
               "/Pathway_and_signatures_results_XtoM.rds"))


# 3- Export counts of pvalues for each set of XtoM----------------
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
fwrite(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets", 
               "/Pvalue_counts_XtoM_SEDAC.csv"))



