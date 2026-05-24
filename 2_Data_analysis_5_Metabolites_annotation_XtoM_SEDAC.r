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
               "apLCMS", 
               "xMSanalyzer",
               "statmod"
               )
# # install xMSanalyzer, and apLCMS, which we will use for metabolites annotation
# # OF NOTE: before installing them, you may need to mannually install ROCS package locally, as well as other dependencies
# # including but are not limited to: xcms, ProtGenerics, sva, mixOmics, WGCNA, GO.db
# # when you finish this, run the remote::install code below finishing installation
# remotes::install_github("kuppal2/xMSanalyzer")
# remotes::install_github("tianwei-yu/apLCMS")

# 1- set work directory----------
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# 2- import untargeted metabolomic mz and retention time information----
fread(str_c(outpath, "/log2.transformed.metabolomics_C18.csv"))  %>% 
    .[, 1:10] -> c18_info
fread(str_c(outpath, "/log2.transformed.metabolomics_HILIC.csv"))  %>% 
    .[, 1:10] -> hilic_info

# import MWAS results
out.TRAPtoM.SEDAC = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
                          "/MWAS_TRAP_Metabolome_SEDAC.rds"))
out.TRAPtoM.mixture = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
                          "/MWAS_TRAP_Metabolome_mixture.rds"))
out.TRAPtoM.mixture %>% list.names()
out.TRAPtoM.mixture$sedac.pm25$c18$coef %>% list.names() -> list.names

groups = c("main", "pm25", "metal", "pm25.full")
f_merge = function(group = "main", lst){
  lst[[str_c("SEDAC_mixture.", group, "_1y")]][["c18"]]    = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$c18$coef[["1y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  lst[[str_c("SEDAC_mixture.", group, "_1y")]][["hilic"]]  = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$hilic$coef[["1y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  lst[[str_c("SEDAC_mixture.", group, "_3y")]][["c18"]]    = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$c18$coef[["3y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  lst[[str_c("SEDAC_mixture.", group, "_3y")]][["hilic"]]  = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$hilic$coef[["3y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  lst[[str_c("SEDAC_mixture.", group, "_5y")]][["c18"]]    = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$c18$coef[["5y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  lst[[str_c("SEDAC_mixture.", group, "_5y")]][["hilic"]]  = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$hilic$coef[["5y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  lst[[str_c("SEDAC_mixture.", group, "_10y")]][["c18"]]   = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$c18$coef[["10y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  lst[[str_c("SEDAC_mixture.", group, "_10y")]][["hilic"]] = out.TRAPtoM.mixture[[str_c("sedac.", group)]]$hilic$coef[["10y"]] %>% dplyr::rename(Components = "Outcome", Pvalue_a = "P.value")
  return(lst)
}

out.TRAPtoM.SEDAC = f_merge("main", lst = out.TRAPtoM.SEDAC)
out.TRAPtoM.SEDAC = f_merge("pm25", lst = out.TRAPtoM.SEDAC)
out.TRAPtoM.SEDAC = f_merge("metal", lst = out.TRAPtoM.SEDAC)
out.TRAPtoM.SEDAC = f_merge("pm25.full", lst = out.TRAPtoM.SEDAC)

# glance
out.TRAPtoM.SEDAC %>% list.names()

# import authentic metabolites information of level 1
fread(str_c(inpath, "/c18 total_confirmed.csv")) -> c18_cfm #neg
fread(str_c(inpath, "/Hilic total_confirmed.csv")) -> hilic_cfm #pos

# FUNCTION to execute annotation for listed files----
F_annotation = function(
  grp_name = "SEDAC_PM2.5_Zinc_1y",
  p.type = "raw", 
  p.threshold = 0.05
  ){
  
  tictoc::tic()

  # assign data
  res = out.TRAPtoM.SEDAC 

  # import data----
  # c18 ----
  res[[grp_name]][["c18"]] %>% 
    dplyr::filter(!is.na(Components)) %>% 
    dplyr::mutate(mid = as.numeric(str_extract(Components, "(?<=V)\\d+")))  -> test1 

  switch(p.type, 
    "raw" = {
      c18_info %>% 
        dplyr::mutate(mid = row_number()) %>% 
        full_join(test1, by = "mid") %>% 
        dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
        dplyr::select(mz, time, Components, Pvalue_a, P_fdr) %>% 
        dplyr::filter(Pvalue_a < p.threshold) -> neg
    }, 
    "fdr" = {
      c18_info %>% 
        dplyr::mutate(mid = row_number()) %>% 
        full_join(test1, by = "mid") %>% 
        dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
        dplyr::select(mz, time, Components, Pvalue_a, P_fdr) %>% 
        dplyr::filter(P_fdr < p.threshold) -> neg
     }
  )# end switch

  # hilic ----
  res[[grp_name]][["hilic"]]  %>% 
    dplyr::filter(!is.na(Components)) %>% 
    dplyr::mutate(mid = as.numeric(str_extract(Components, "(?<=V)\\d+"))) -> test2

  switch(p.type, 
    "raw" = {
      hilic_info %>% 
        dplyr::mutate(mid = row_number()) %>% 
        full_join(test2, by = "mid") %>% 
        dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
        dplyr::select(mz, time, Components, Pvalue_a, P_fdr) %>% 
        dplyr::filter(Pvalue_a < p.threshold) -> pos
    }, 
    "fdr" = {
      hilic_info %>% 
        dplyr::mutate(mid = row_number()) %>% 
        full_join(test2, by = "mid") %>% 
        dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
        dplyr::select(mz, time, Components, Pvalue_a, P_fdr) %>% 
        dplyr::filter(P_fdr < p.threshold) -> pos
    }
  )# end switch
  
  # negative mode (c18)----
  ovlp.neg = xMSanalyzer::find.Overlapping.mzs(
              neg[ ,1:2], 
              map_dfc(c18_cfm[, 5:6] , ~as.numeric(.x)),
              mz.thresh = 10, 
              # time.thresh = 60, we usually don't use this, since rt would be influenced by various external conditions
              alignment.tool = 'apLCMS')
  if(nrow(ovlp.neg) == 0) {
    print(str_c("NO authentic metabolites have been annotated for positive mode in ", grp_name))
    neg.ant = data.table()
  } else {
    neg[ovlp.neg$index.A,] %>% 
      cbind(
      c18_cfm[ovlp.neg$index.B, c(5, 6, 3, 4, 2)] %>% 
        set_names(c("mz_confirm", "rt_confirm", "HMDB_monoisotopic_mass", "Adduct", "Metabolite_name"))
      ) %>% 
      dplyr::mutate(Source = str_extract(grp_name, "^(Rline|SEDAC)"), #Rline or SEDAC
                    Exposure = str_extract(grp_name, "(?<=(Rline_|SEDAC_))[a-zA-Z0-9._]+(?=_\\d+y$)"), # air pollution components
                    Window = str_extract(grp_name, "(?<=_)[0-9]+y$"), # time windows
                    Mode = "Negative") -> neg.ant
  }
    
  # positive mode (hilic)----
  ovlp.pos = xMSanalyzer::find.Overlapping.mzs(
              pos[ ,1:2], 
              map_dfc(hilic_cfm[, 5:6] , ~as.numeric(.x)),
              mz.thresh = 10, 
              # time.thresh = 60, we usually don't use this, since rt would be influenced by various external conditions
              alignment.tool = 'apLCMS')
  if(nrow(ovlp.pos) == 0) {
    print(str_c("NO authentic metabolites have been annotated for positive mode in ", grp_name))
    pos.ant = data.table()
  } else {
    pos[ovlp.pos$index.A,] %>% 
      cbind(
      hilic_cfm[ovlp.pos$index.B, c(5, 6, 3, 4, 2)] %>% 
        set_names(c("mz_confirm", "rt_confirm", "HMDB_monoisotopic_mass", "Adduct", "Metabolite_name"))
      ) %>% 
      dplyr::mutate(Source = str_extract(grp_name, "^(Rline|SEDAC)"), #Rline or SEDAC
                    Exposure = str_extract(grp_name, "(?<=(Rline_|SEDAC_))[a-zA-Z0-9._]+(?=_\\d+y$)"), # air pollution components
                    Window = str_extract(grp_name, "(?<=_)[0-9]+y$"), # time windows
                    Mode = "Positive") -> pos.ant
  } 
  # combine negative and positive results together
  rbind(neg.ant, pos.ant) -> ant

  return(ant)
}

# 4 Execute function----
# Create output folder
if(!file.exists(str_c(outpath,"/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets"))){
  dir.create(str_c(outpath,"/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets"))
}

## 4.1 execute: raw p <0.05----
# SEDAC
sedac_raw0.05 = purrr::pmap(list(
                      out.TRAPtoM.SEDAC %>% list.names(), 
                      rep("raw", times = length(out.TRAPtoM.SEDAC %>% list.names())), 
                      rep(0.05, times = length(out.TRAPtoM.SEDAC %>% list.names()))),
                      F_annotation)  %>% rbindlist
sedac_raw0.05 %>% 
  fwrite(str_c(outpath,"/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
              "/Annotated_metabolites_XtoM_SEDAC_thresholds_mz10ppm_rawp0.05.csv"))
sedac_raw0.05$Metabolite_name %>% unique() %>% length()

sedac_raw0.05 %>% 
  dplyr::select(mz, time, Mode) %>% unique() %>% nrow()

## 4.2 execute: fdr p <0.2----
# SEDAC
sedac_fdr0.2 = purrr::pmap(list(
                      out.TRAPtoM.SEDAC %>% list.names(), 
                      rep("fdr", times = length(out.TRAPtoM.SEDAC %>% list.names())), 
                      rep(0.2, times = length(out.TRAPtoM.SEDAC %>% list.names()))),
                      F_annotation)  %>% rbindlist
sedac_fdr0.2 %>% 
  fwrite(str_c(outpath,"/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
              "/Annotated_metabolites_XtoM_SEDAC_thresholds_mz10ppm_fdr0.2.csv"))
sedac_fdr0.2$Metabolite_name %>% unique() %>% length()


# extract list
sedac.exposures = out.TRAPtoM.SEDAC %>% list.names() %>% str_extract("(?<=SEDAC_)[a-zA-Z0-9._]+(?=_\\d+y$)") %>% unique()

# data = sedac_fdr0.2
# p.class = "fdr0.2"
# 6- Draw plots for annotated chemicals using diffrent p-thresholds----
f_sedac = function(data, p.class){
  tibble(Exposure = sedac.exposures, N = rep(0,length(sedac.exposures))) %>% 
    full_join(tibble(P_threshold = p.class, N = 0), by = "N", relationship = "many-to-many") %>% 
    full_join(tibble(Window = c("1y", "3y", "5y", "10y"), N = rep(0,4)), by = "N", relationship = "many-to-many") %>% 
    dplyr::select(Exposure, Window, N, P_threshold) -> tmp
  if(nrow(data)!=0){
    data %>% 
      group_by(Exposure, Window) %>% 
      summarise(N = unique(Metabolite_name) %>% length())  %>% 
      dplyr::mutate(P_threshold = p.class) -> tmp2
    tmp %>% 
      rows_update(tmp2, by = c("Exposure", "Window"))
  }else if (nrow(data)==0) {
    tmp
  }
}
# execute
purrr::map2(list(sedac_raw0.05, sedac_fdr0.2), 
            list("raw0.05", "fdr0.2"),
            f_sedac) %>% rbindlist() -> plotdata

plotdata %<>% 
  dplyr::mutate(P_threshold = factor(P_threshold, levels = c("raw0.05", "fdr0.2")), 
                Exposure = factor(Exposure, levels = sedac.exposures), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y"))) %>% 
  arrange(P_threshold, Exposure, Window)

## 6.1-plot ----
plotdata$Exposure %>% unique() -> plt.list

plotdata %>% 
  dplyr::filter(Exposure  %in% plt.list[c(1,2,3,19)]) %>% 
  ggplot(aes(P_threshold, Window)) +
  geom_tile(aes(fill = N)) + 
  scale_fill_gradient(low = "white", high = "blue") + 
  facet_grid(Exposure~.) + 
  theme(axis.ticks.y=element_blank(),
     axis.title.y = element_text(vjust = +3),
     axis.text.y = element_text(size=16), 
     legend.position = "bottom") +
  ggtitle("Annotated number of metabolites in different p thresholds that associated with TRAP")
ggsave(str_c(outpath, "/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets/Annotated_metabolites_counts_SEDAC_main", ".png"), height = 12, width = 10, dpi = 600)

plotdata %>% 
  dplyr::filter(!(Exposure  %in% plt.list[c(1,2,3,19)])) %>%
  ggplot(aes(P_threshold, Window)) +
  geom_tile(aes(fill = N)) + 
  scale_fill_gradient(low = "white", high = "blue") + 
  facet_grid(Exposure~.) + 
  theme(axis.ticks.y=element_blank(),
     axis.title.y = element_text(vjust = +3),
     axis.text.y = element_text(size=16), 
     legend.position = "bottom") +
  ggtitle("Annotated number of metabolites in different p thresholds that associated with TRAP")
ggsave(str_c(outpath, "/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets/Annotated_metabolites_counts_SEDAC_pm25", ".png"), height = 30, width = 10, dpi = 600)

# export table
plotdata %>% 
  pivot_wider(names_from = "P_threshold", values_from = "N") %>% 
  fwrite(str_c(outpath,"/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
              "/Annotated_metabolites_counts_SEDAC.csv"))


# plot fdr q<0.2-----------
sedac_fdr0.2 = 
  fread(str_c(outpath,"/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
              "/Annotated_metabolites_XtoM_SEDAC_thresholds_mz10ppm_fdr0.2.csv"))
sedac_fdr0.2$Metabolite_name %>% unique() %>% length()

sedac_fdr0.2 %>% 
  dplyr::filter(!(Exposure  %in% c("O3", "NO2", "mixture.main"))) -> plt
lst = plt$Exposure %>% unique()

plt %<>% 
  dplyr::mutate(Exposure = factor(Exposure, levels = lst[c(9,8,7,1,2,3,4,5,6)]), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y"))) %>% 
  arrange(Exposure, Window)

plt %>% 
  group_by(Exposure, Window) %>% 
      summarise(N = unique(Metabolite_name) %>% length()) %>% 
  ungroup() -> plt2
  
plt2  %>% 
  ggplot(aes(Exposure, Window)) +
  geom_tile(aes(fill = N)) + 
  scale_fill_gradient(low = "white", high = "blue", na.value = "grey") + 
  theme_bw() + 
  theme(axis.ticks.y=element_blank(),
     axis.title.y = element_text(vjust = +3),
     axis.text.y = element_text(size=16), 
     legend.position = "bottom")
ggsave(str_c(outpath, "/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets/Annotated_metabolites_counts_SEDAC_pm25(fdr q0.2)", ".tiff"), height = 10, width = 10, dpi = 600)





######################
###export files for mzmine-----------

# fdr q<0.2
sedac_fdr0.2 =
  fread(str_c(outpath,"/DataAnalysis_3_Metabolite_annotation_XtoM/Sensitive_analyses_using maximum datasets", 
              "/Annotated_metabolites_XtoM_SEDAC_thresholds_mz10ppm_fdr0.2.csv"))
sedac_fdr0.2 %>% 
  rbind(rline_fdr0.2) %>% 
  dplyr::select(mz, time, Components, Mode, Metabolite_name) %>% 
  unique() %>% 
  group_by(mz, time, Components, Mode) %>% 
  dplyr::mutate(n = 1:n()) %>% 
  pivot_wider(names_from = "n", values_from = "Metabolite_name", names_prefix = "metabolite_") %>% 
  fwrite(str_c(outpath, "/0_mzmine/Annotated_metabolites_fdr0.2_new.csv"))

# -----------------------------------
# find metabolites that exist in new dataset and not included in previous dataset 
old = fread(str_c(outpath, "/0_mzmine/Annotated_metabolites_fdr0.2.csv"))
new = fread(str_c(outpath, "/0_mzmine/Annotated_metabolites_fdr0.2_new.csv"))
old %>% 
  dplyr::select(mz, time, Components, Mode) %>% 
  unique() %>% 
  dplyr::mutate(flag = "old")  -> tmp
new %>% 
  left_join(tmp, by = c("mz", "time", "Components", "Mode")) %>% 
  dplyr::filter(is.na(flag)) -> new2
# export new metabolites that have never been EIC-extracted
new2 %>% 
  fwrite(str_c(outpath, "/0_mzmine/Annotated_metabolites_fdr0.2_update_set.csv"))
# ----------------------------------
