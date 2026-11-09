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

##############################################################
# Describe confirmed and validated metabolite levels across AD stages
######################################################################

# In hypothesis-generating cohort------------------------------------
# set work directory
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# load data
load(str_c(outpath, "/!Data_for_Data_Analysis_1.rdata"))
ls()

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

# get the dataset to be used for descriptive analysis 
df.c18 %>% 
  dplyr::filter(!is.na(df.c18[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1)) -> tmp.c18

df.hilic %>% 
  dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1)) -> tmp.hilic

# read metabolite list to be used
metabolites <- fread(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_fdr0.2_final.csv")) %>% 
              dplyr::select(mz, time, Components, Mode, metabolite_1)

# restricted to those replicated
metabolites %>% 
  dplyr::filter(str_detect(metabolite_1, "COA|METHIONINE")) -> replicated


# loop function

var = "V7277"
mode = "negative"

BoxPlot = function(var, mode){
  switch(mode,
         "negative" = {
          tmp <- tmp.c18
          },  
         "positive" = {
          tmp <- tmp.hilic
          }
        )
  
  tmp %>% 
    dplyr::select(all_of(c("ABC_cat", "Braak_cat", "CERAD_cat")), # AD pathology markers
                  all_of(var))  %>% 
    pivot_longer(cols = all_of(c("ABC_cat", "Braak_cat", "CERAD_cat")), names_to = "ADmarker", values_to = "value") %>% 
    dplyr::mutate(ADmarker = factor(ADmarker, 
                                    levels = tmp %>% 
                                    dplyr::select(all_of(c("ABC_cat", "Braak_cat", "CERAD_cat"))) %>% 
                                    colnames()))  -> plotdata
  # extract maximum range of plot
  max_y_range = max(plotdata[[var]]) - min(plotdata[[var]])

  plotdata %>% 
    dplyr::select(any_of(var), ADmarker, value)  %>% 
    set_names(c("var", "ADmarker", "value")) %>% 
    dplyr::filter(!is.na(var)) -> a
  
  a  %<>% 
    dplyr::mutate(value = case_when(value %in% c("High", "Higher stage", "Frequent") ~ "1",
                                    value %in% c("Low", "Lower stage", "Sparse") ~ "0",
                                    TRUE ~ NA)) %>% 
    dplyr::mutate(value = as.integer(value))

  # stats
  q1 = quantile(plotdata[[var]], probs = c(0.25))
  q3 = quantile(plotdata[[var]], probs = c(0.75))
  iqr = q3 - q1

  r1 = q1 - 3 * iqr
  r2 = q3 + 3 * iqr
  # 

  a %>% 
    dplyr::mutate(value = factor(value, levels = c("0", "1"))) %>%
    ggplot(aes(x = value, y = var, fill = value, group = value)) + 
     geom_half_violin(aes(fill=value, color=value),trim=FALSE, side="r") + # half violin
      # geom_half_dotplot(aes(x=as.numeric(value) - 0.07, fill=value, color=NA),
      #                   alpha = 0.2, method="histodot", 
      #                   #binwidth = 0.02,
      #                   dotsize =  0.3*(max(plotdata[[var]]) - min(plotdata[[var]]))/max_y_range,
      #                   stackdir="down", stackratio = 0.28, 
      #                   position = position_jitter(height = .03, width = .005)) + #half dot plot
      geom_boxplot(color="black",fill="white", outlier.colour = NA, width=0.1) + # box
      stat_compare_means(comparisons = list(c("0", "1")), 
                         method = "wilcox.test", 
                         method.args = list(alternative = "two.sided"),
                         label.y = r2 * 0.98,
                         tip.length = 0.01,  
                         label = "p.signif") + 
      guides(fill=FALSE, color=FALSE) + 
      theme( #plot parameters
          panel.border = element_rect(fill=NA,color="black", size=0.7, linetype="solid"),
          panel.background=element_blank(),
          axis.text = element_text(size = 14,  color = "black", face = "plain"),
          axis.title = element_text(size = 15, color = "black", face = "plain"), 
          strip.text = element_text(size = 14, color = "black", face = "plain")) +
      coord_cartesian() +
    scale_fill_manual(values=c("#7cadca", "#b1564c")) + 
    scale_color_manual(values=c("#7cadca", "#b1564c")) + 
    facet_grid(~ADmarker ,scales = "free_x") +  
      labs(
        title = str_c("Distributions of ", var, " among AD pathology marker groups"),
        y = "log2-transformed metabolite levels",
        x = glue("AD pathology stages (0-early, 1-advanced)")) + 
    ylim(r1, r2*1.02)

  ggsave(str_c(outpath, "/Visualizations/0_SM_figs", 
         "/Extended_Fig_metabolite levels_[", var, "_", mode, "]_among AD markers.png"), 
         width = 16, 
         height = 10,
         dpi = 600)

  print(str_c("Distributions of metabolite: ", var, "_", mode, " accomplished!"))
}

# execute
map2(replicated$Components, replicated$Mode, BoxPlot)


var = "V7277"
mode = "negative"

# describle relations with PMI, using scatter plots
DotPlot = function(var, mode){
  switch(mode,
         "negative" = {
          tmp <- tmp.c18
          },  
         "positive" = {
          tmp <- tmp.hilic
          }
        )
  
  tmp %>% 
    dplyr::select(all_of(c("PMI_MY")), # PMI
                  all_of(var)) %>% 
    dplyr::rename(var = var) -> plotdata
  
  # extract maximum range of plot
  max_y_range = max(plotdata$var) - min(plotdata$var)

  plotdata %>%
    ggplot(aes(x = PMI_MY, y = var)) +
    geom_point(color = "#7cadca", size = 2.5, alpha = 0.8) +                  # points
    geom_smooth(method = "lm", se = TRUE, color = "#b1564c",
                fill = "#b1564c", alpha = 0.15, linewidth = 0.9) +            # fit + 95% CI
    stat_cor(method = "spearman",                                             # or "pearson"
             cor.coef.name = "rho",
             label.x.npc = "right", label.y.npc = "bottom", 
             hjust = 1,
             vjust = 0,
             size = 5) +
    theme( # plot parameters
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7, linetype = "solid"),
      panel.background = element_blank(),
      axis.text = element_text(size = 14, color = "black", face = "plain"),
      axis.title = element_text(size = 15, color = "black", face = "plain"),
      strip.text = element_text(size = 14, color = "black", face = "plain")) +
    labs(
      title = str_c("Association between ", var, " and PMI"),
      y = "log2-transformed metabolite levels",
      x = "Post-mortem interval (PMI)")

  ggsave(str_c(outpath, "/Visualizations/0_SM_figs", 
         "/Extended_Fig_metabolite levels_[", var, "_", mode, "]_with PMI.png"), 
         width = 10, 
         height = 10,
         dpi = 600)

  print(str_c("Distributions of metabolite: ", var, "_", mode, " accomplished!"))
}

# execute
map2(replicated$Components, replicated$Mode, DotPlot)




############################################
rm(list=(ls()))
## describe same things in replication dataset
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs/1_Internal_validation")

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
fread(str_c(outpath, "/1_Met2025_metabolites_to_be_modelled.csv")) %>%
  dplyr::select(-V1) %>% 
  as_tibble() %>%
  dplyr::mutate(V1 = row_number()) -> m.data

dim(m.data)

m.data %>% 
  dplyr::select(14:ncol(m.data))  %>% 
  t() %>% 
  as_tibble() %>% 
  dplyr::mutate(Sample_ID = colnames(m.data[, 14:ncol(m.data)]), 
                GUID = case_when(str_detect(Sample_ID, "^ID_") ~ str_replace_all(Sample_ID, "ID_", ""), 
                                TRUE ~ Sample_ID)) %>% 
  dplyr::select(GUID, everything(), -Sample_ID) %>% 
  dplyr::filter(!str_detect(GUID, "^q4|^nist")) -> meta.data 

dim(meta.data)
# impute missings as mim/2
meta.data %<>%
  dplyr::mutate(across(all_of(matches("V\\d+")), ~ifelse(is.na(.x), min(.x/2, na.rm = TRUE), .x)))

meta.data$GUID = as.integer(meta.data$GUID)

# merge data
demo %>% 
  left_join(ap.data, by = "Case_Number") -> tmp

tmp %>% 
  dplyr::rename(GUID = "GUID.x") %>%
  left_join(meta.data, by = "GUID") %>% 
  as_tibble() -> df.vali

# 
df.vali %>% 
  dplyr::filter(!is.na(V1)) %>% # those have metabolomics data -- 347 in total
  dplyr::filter(!is.na(.[["SEDAC_PM25NO3_1y"]])) -> tmp # those have air pollution data -- 300 in total

# import those have been analyzed in discovery phase
df.dis = fread(str_c(work.dir, "/3_Outputs/Cleaned_TRAP_data_for_brain_donor_N162.csv"))
# GUID in discovery cohort
df.dis %>% 
  dplyr::filter(!is.na(SEDAC_PM25NO3_1y)) %>% 
  .$GUID -> guid.dis

# final data for validation
tmp %>% 
  dplyr::filter(!(GUID %in% guid.dis)) -> df.validation #143 obs for validation

colnames(df.validation)

df.validation %<>% 
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

m.data %>% 
  dplyr::mutate(new_component = str_c("V", row_number())) %>% 
  dplyr::select(new_component, mode:time) %>% 
  dplyr::filter(str_detect(metabolite, "CoA|Methionine")) -> replicated

var = "V7277"
mode = "Negative"

BoxPlot_rep = function(var, mode){

  tmp <- df.validation
  var2 <- replicated$new_component[which(replicated$Components == var)]

  tmp %>% 
    dplyr::select(all_of(c("ABC_cat", "Braak_cat", "CERAD_cat")), # AD pathology markers
                  all_of(var2))  %>% 
    pivot_longer(cols = all_of(c("ABC_cat", "Braak_cat", "CERAD_cat")), names_to = "ADmarker", values_to = "value") %>% 
    dplyr::mutate(ADmarker = factor(ADmarker, 
                                    levels = tmp %>% 
                                    dplyr::select(all_of(c("ABC_cat", "Braak_cat", "CERAD_cat"))) %>% 
                                    colnames())) %>% 
    dplyr::filter(!is.na(value))  -> plotdata
  # extract maximum range of plot
  max_y_range = max(plotdata[[var2]]) - min(plotdata[[var2]])

  plotdata %>% 
    dplyr::select(any_of(var2), ADmarker, value)  %>% 
    set_names(c("var2", "ADmarker", "value")) %>% 
    dplyr::filter(!is.na(var2)) -> a
  
  a  %<>% 
    dplyr::mutate(value = case_when(value %in% c("High", "Higher stage", "Frequent") ~ "1",
                                    value %in% c("Low", "Lower stage", "Sparse") ~ "0",
                                    TRUE ~ NA)) %>% 
    dplyr::mutate(value = as.integer(value))

  # stats
  q1 = quantile(plotdata[[var2]], probs = c(0.25))
  q3 = quantile(plotdata[[var2]], probs = c(0.75))
  iqr = q3 - q1

  r1 = q1 - 3 * iqr
  r2 = q3 + 3 * iqr
  # 

  a %>% 
    dplyr::mutate(value = factor(value, levels = c("0", "1"))) %>%
    ggplot(aes(x = value, y = var2, fill = value, group = value)) + 
     geom_half_violin(aes(fill=value, color=value),trim=FALSE, side="r") + # half violin
      # geom_half_dotplot(aes(x=as.numeric(value) - 0.07, fill=value, color=NA),
      #                   alpha = 0.2, method="histodot", 
      #                   #binwidth = 0.02,
      #                   dotsize =  0.3*(max(plotdata[[var]]) - min(plotdata[[var]]))/max_y_range,
      #                   stackdir="down", stackratio = 0.28, 
      #                   position = position_jitter(height = .03, width = .005)) + #half dot plot
      geom_boxplot(color="black",fill="white", outlier.colour = NA, width=0.1) + # box
      stat_compare_means(comparisons = list(c("0", "1")), 
                         method = "wilcox.test", 
                         method.args = list(alternative = "two.sided"),
                         label.y = r2 * 0.98,
                         tip.length = 0.01,  
                         label = "p.signif") + 
      guides(fill=FALSE, color=FALSE) + 
      theme( #plot parameters
          panel.border = element_rect(fill=NA,color="black", size=0.7, linetype="solid"),
          panel.background=element_blank(),
          axis.text = element_text(size = 14,  color = "black", face = "plain"),
          axis.title = element_text(size = 15, color = "black", face = "plain"), 
          strip.text = element_text(size = 14, color = "black", face = "plain")) +
      coord_cartesian() +
    scale_fill_manual(values=c("#7cadca", "#b1564c")) + 
    scale_color_manual(values=c("#7cadca", "#b1564c")) + 
    facet_grid(~ADmarker ,scales = "free_x") +  
      labs(
        title = str_c("Distributions of ", var, " among AD pathology marker groups(replication data)"),
        y = "log2-transformed metabolite levels",
        x = glue("AD pathology stages (0-early, 1-advanced)")) + 
    ylim(r1, r2*1.02)

  ggsave(str_c(work.dir, "/3_Outputs/Visualizations/0_SM_figs", 
         "/Extended_Fig_metabolite levels_[", var, "_", mode, "]_among AD markers(replication data).png"), 
         width = 16, 
         height = 10,
         dpi = 600)

  print(str_c("Distributions of metabolite: ", var, "_", mode, " accomplished!"))
}

# execute
map2(replicated$Components, replicated$mode, BoxPlot_rep)


# describle relations with PMI, using scatter plots
DotPlot_rep = function(var, mode){

  tmp <- df.validation
  var2 <- replicated$new_component[which(replicated$Components == var)]

  tmp %>% 
    dplyr::select(all_of(c("PMI_MY")),
                  all_of(var2)) %>% 
    dplyr::rename(var2 = var2)  -> plotdata

  plotdata %>%
    ggplot(aes(x = PMI_MY, y = var2)) +
    geom_point(color = "#7cadca", size = 2.5, alpha = 0.8) +                  # points
    geom_smooth(method = "lm", se = TRUE, color = "#b1564c",
                fill = "#b1564c", alpha = 0.15, linewidth = 0.9) +            # fit + 95% CI
    stat_cor(method = "spearman",                                             # or "pearson"
             cor.coef.name = "rho",
             label.x.npc = "right", label.y.npc = "bottom", 
             hjust = 1,
             vjust = 0,
             size = 5) +
    theme( # plot parameters
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7, linetype = "solid"),
      panel.background = element_blank(),
      axis.text = element_text(size = 14, color = "black", face = "plain"),
      axis.title = element_text(size = 15, color = "black", face = "plain"),
      strip.text = element_text(size = 14, color = "black", face = "plain")) +
    labs(
      title = str_c("Association between ", var, " and PMI"),
      y = "log2-transformed metabolite levels",
      x = "Post-mortem interval (PMI)")

  ggsave(str_c(work.dir, "/3_Outputs/Visualizations/0_SM_figs", 
         "/Extended_Fig_metabolite levels_[", var, "_", mode, "]_with PMI(replication data).png"), 
         width = 10, 
         height = 10,
         dpi = 600)

  print(str_c("Distributions of metabolite: ", var, "_", mode, " accomplished!"))
}


# execute
map2(replicated$Components, replicated$mode, DotPlot_rep)





df.validation$ABC_new %>% table(useNA = "always")
df.validation$Braak_Stage_new %>% table(useNA = "always")
df.validation$CERAD_new %>% table(useNA = "always")


######################################################################
######################################################################
######################################################################
######################################################################

# describe metabolites levels among subgroups
# (to be used as supplementary resulgts along with stratified analysis)

# In hypothesis-generating cohort------------------------------------
rm(list=(ls()))
# set work directory
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# load data
load(str_c(outpath, "/!Data_for_Data_Analysis_1.rdata"))
ls()

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
                                      TRUE ~ NA)) %>% 
  dplyr::mutate(Age_cat = factor(Age_cat, levels = c("<70 yr", ">=70 yr")),
                sex = factor(sex, levels = c("Male", "Female")),
                ApoE_cat = factor(ApoE_cat, levels = c("E4 absence", "E4 existance")),
                ABC_cat = factor(ABC_cat, levels = c("Low", "High")), 
                Braak_cat = factor(Braak_cat, levels = c("Lower stage", "Higher stage")), 
                CERAD_cat = factor(CERAD_cat, levels = c("Sparse", "Frequent")))

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
                                      TRUE ~ NA)) %>% 
  dplyr::mutate(Age_cat = factor(Age_cat, levels = c("<70 yr", ">=70 yr")),
                sex = factor(sex, levels = c("Male", "Female")),
                ApoE_cat = factor(ApoE_cat, levels = c("E4 absence", "E4 existance")),
                ABC_cat = factor(ABC_cat, levels = c("Low", "High")), 
                Braak_cat = factor(Braak_cat, levels = c("Lower stage", "Higher stage")), 
                CERAD_cat = factor(CERAD_cat, levels = c("Sparse", "Frequent")))

# get the dataset to be used for descriptive analysis 
df.c18 %>% 
  dplyr::filter(!is.na(df.c18[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1)) -> tmp.c18

df.hilic %>% 
  dplyr::filter(!is.na(df.hilic[[str_c("SEDAC_PM25NH4_1y")]]) & !is.na(V1)) -> tmp.hilic

# read metabolite list to be used
metabolites <- fread(str_c(outpath, "/0_mzmine/mzmine_fdr0.2/Youran_Mengyuan_Donghai", 
                                      "/0_Confirmed_metabolites_fdr0.2_final.csv")) %>% 
              dplyr::select(mz, time, Components, Mode, metabolite_1)


## describe exposure levels among subgroups ----------------------
tmp.c18 %>% 
  dplyr::select(contains(c("_O3_", "PM2.5_Lead", "PM2.5_Silicon", "PM2.5_Vanadium", "PM25SO4"))) %>% 
  dplyr::select(!contains("_15y")) %>% 
  colnames() -> expo.list


# loop function 
var = "ABC_cat"
mode = "negative" # it doesnt matter since we are comparing exposure levels among subgroups

BoxPlot = function(var, mode){
  switch(mode,
         "negative" = {
          tmp <- tmp.c18
          },  
         "positive" = {
          tmp <- tmp.hilic
          }
        )
  
  tmp %>% 
    dplyr::select(all_of(expo.list), # air pollutants
                  Age_cat, sex, ApoE_cat, ABC_cat, Braak_cat, CERAD_cat) %>% 
    pivot_longer(cols = all_of(expo.list), names_to = "AirPollutant", values_to = "value") %>% 
    dplyr::mutate(AirPollutant = factor(AirPollutant, 
                               levels = tmp %>% 
                               dplyr::select(all_of(expo.list)) %>% 
                               colnames()))  -> plotdata

  # extract maximum range of plot
  plotdata %>% 
    dplyr::mutate(subgrp = str_extract(AirPollutant, "(?<=SEDAC_)[0-9a-zA-Z\\.\\_]+(?=_)"), 
                  tms = str_extract(AirPollutant, "(?<=_)\\d+y$")) ->  a
  max_y_range = max(a$value, na.rm = TRUE) - min(a$value, na.rm = TRUE)

  a %>% 
    dplyr::select(any_of(var), AirPollutant, value)  %>% 
    set_names(c("var", "AirPollutant", "value")) %>% 
    dplyr::filter(!is.na(var)) -> a
  
  # stats
  q1 = quantile(plotdata[["value"]], probs = c(0.25), na.rm = TRUE)
  q3 = quantile(plotdata[["value"]], probs = c(0.75), na.rm = TRUE)
  iqr = q3 - q1

  r1 = q1 - 1.5 * iqr
  r2 = q3 + 1.5 * iqr

  comparisons = combn(unique(a$var), 2, simplify = FALSE)

  a %>% 
    dplyr::mutate(subgrp = str_extract(AirPollutant, "(?<=SEDAC_)[0-9a-zA-Z\\.\\_]+(?=_)"), 
                  tms = str_extract(AirPollutant, "(?<=_)\\d+y$"), 
                  subgrp = str_replace(subgrp, "PM2.5_|PM25", "")) %>%
    dplyr::mutate(   
      tms    = factor(tms,    levels = c("1y", "3y", "5y", "10y")),
      subgrp = factor(subgrp, levels = c("O3", "SO4", "Lead", "Silicon", "Vanadium"))
    ) %>% 
    ggplot(aes(x = var, y = log2(value), fill = var, group = var)) + 
     geom_half_violin(aes(fill=var, color=var),trim=FALSE, side="r") + # half violin
      geom_boxplot(color="black",fill="white", outlier.colour = NA, width=0.1) + # box
      stat_compare_means(comparisons = list(levels(a$var)), 
                         label = "p.signif") + 
      guides(fill=FALSE, color=FALSE) + 
      theme( #plot parameters
          panel.border = element_rect(fill=NA, color="black", size=0.7, linetype="solid"),
          panel.background=element_blank(),
          axis.text = element_text(size = 14,  color = "black", face = "plain"),
          axis.title = element_text(size = 15, color = "black", face = "plain"), 
          strip.text = element_text(size = 14, color = "black", face = "plain")) +
      coord_cartesian() +
    scale_fill_manual(values=c("#7cadca", "#b1564c", "#E69F00")) + 
    scale_color_manual(values=c("#7cadca", "#b1564c", "#E69F00")) + 
    facet_grid(subgrp~tms ,scales = "free_y") +  
      labs(
        title = str_c("Distributions of air pollutants of various exposure windows among group ", var),
        y = "log2-transformed exposure concentrations",
        x = var)
  
  ggsave(str_c(outpath, "/Visualizations/0_SM_figs", 
         "/Extended_Fig_air pollution level_comparison_among_[", var, "]", ".png"), 
         width = 16, 
         height = 16,
         dpi = 600)

  print(str_c("Distributions of air pollution of various exposure windows among group: ", var, " accomplished!"))

}


# execute 
map2(c("Age_cat", "sex", "ApoE_cat", "ABC_cat", "Braak_cat", "CERAD_cat"),
     rep("negative", 6), 
     BoxPlot)


## describe metabolite levels among subgroups ----------------------
# read metabolite list to be used
metabolites <- readxl::read_xlsx(str_c(work.dir, "/0_Manuscript/tables", 
                                      "/PN_Identified metabolites (Level-1 evidence).xlsx")) %>% 
              dplyr::select(all_of(c("m/z", "rt (second)", "Components", "Mode", "Brief name"))) %>% 
              set_names(c("mz", "time", "Components", "Mode", "metabolite_name"))

# loop function 
var = "ABC_cat"

BoxPlot = function(var){

  meta.list.c18 <- metabolites %>% 
    dplyr::filter(Mode == "Negative") %>% 
    .$Components
  meta.list.hilic <- metabolites %>% 
    dplyr::filter(Mode == "Positive") %>% 
    .$Components

  tmp <- tmp.c18 %>% 
          dplyr::select(all_of(meta.list.c18), # c18 features
                  Age_cat, sex, ApoE_cat, ABC_cat, Braak_cat, CERAD_cat) %>% 
          bind_cols(tmp.hilic %>% 
                    dplyr::select(all_of(meta.list.hilic))
                    )

  tmp %>% 
    pivot_longer(cols = all_of(c(meta.list.c18, meta.list.hilic)), names_to = "Metabolites", values_to = "value") %>%
    left_join(metabolites %>% dplyr::select(Components, Mode, metabolite_name), by = c("Metabolites" = "Components")) %>%
    dplyr::mutate(column = case_when(Mode == "Negative" ~ "c18",
                                     Mode == "Positive" ~ "hilic"),
                  Metabolites = glue("{metabolite_name} ({column})")) %>%
    arrange(Metabolites) %>% 
    dplyr::mutate(Metabolites = factor(Metabolites, 
                               levels = .$Metabolites %>% unique()))  -> plotdata

  max_y_range = max(plotdata$value, na.rm = TRUE) - min(plotdata$value, na.rm = TRUE)

  plotdata %>% 
    dplyr::select(any_of(var), Metabolites, value)  %>% 
    set_names(c("var", "Metabolites", "value")) %>% 
    dplyr::filter(!is.na(var)) -> a
  
  a_trim <- a %>%
    group_by(Metabolites, var) %>%            # per panel AND per box
    mutate(
      q1  = quantile(value, 0.25, na.rm = TRUE),
      q3  = quantile(value, 0.75, na.rm = TRUE),
      iqr = q3 - q1
    ) %>%
    filter(value >= q1 - 3 * iqr,
           value <= q3 + 3 * iqr) %>%
    ungroup()
  
  a_trim %>% 
    ggplot(aes(x = var, y = value, fill = var, group = var)) + 
     geom_half_violin(aes(fill=var, color=var),trim=FALSE, side="r") + # half violin
      geom_boxplot(color="black",fill="white", outlier.colour = NA, width=0.1) + # box
      stat_compare_means(comparisons = list(levels(a$var)), 
                         label = "p.signif") + 
      guides(fill=FALSE, color=FALSE) + 
      theme( #plot parameters
          panel.border = element_rect(fill=NA, color="black", size=0.7, linetype="solid"),
          panel.background=element_blank(),
          axis.text = element_text(size = 14,  color = "black", face = "plain"),
          axis.title = element_text(size = 15, color = "black", face = "plain"), 
          strip.text = element_text(size = 14, color = "black", face = "plain")) +
      coord_cartesian() +
    scale_fill_manual(values=c("#7cadca", "#b1564c", "#E69F00")) + 
    scale_color_manual(values=c("#7cadca", "#b1564c", "#E69F00")) + 
    facet_wrap(~Metabolites , scales = "free_y", nrow = 5) +  
      labs(
        title = str_c("Distributions of metabolite levels among group ", var),
        y = "log2-transformed metabolite levels",
        x = var)
  
  ggsave(str_c(outpath, "/Visualizations/0_SM_figs", 
         "/Extended_Fig_metabolite level_comparison_among_[", var, "]", ".png"), 
         width = 16 ,
         height = 16,
         dpi = 600)

  print(str_c("Distributions of metabolite levels among group: ", var, " accomplished!"))

}


# execute 
purrr::map(c("Age_cat", "sex", "ApoE_cat", "ABC_cat", "Braak_cat", "CERAD_cat"),
     BoxPlot)



# ----------------------------------------------------

# loop function 
var = "PMI_MY"

DotPlot = function(var){

  meta.list.c18 <- metabolites %>% 
    dplyr::filter(Mode == "Negative") %>% 
    .$Components
  meta.list.hilic <- metabolites %>% 
    dplyr::filter(Mode == "Positive") %>% 
    .$Components

  tmp <- tmp.c18 %>% 
          dplyr::select(all_of(meta.list.c18), # c18 features
                  PMI_MY) %>% 
          bind_cols(tmp.hilic %>% 
                    dplyr::select(all_of(meta.list.hilic))
                    ) %>% 
          dplyr::select(PMI_MY, everything())

  tmp %<>% 
    dplyr::mutate(across(V646:V5989, ~{
      min = min(., na.rm = TRUE)
      ifelse(.== min, NA, .)
    }))

  tmp %>% 
    pivot_longer(cols = all_of(c(meta.list.c18, meta.list.hilic)), names_to = "Metabolites", values_to = "value") %>%
    left_join(metabolites %>% dplyr::select(Components, Mode, metabolite_name), by = c("Metabolites" = "Components")) %>%
    dplyr::mutate(column = case_when(Mode == "Negative" ~ "c18",
                                     Mode == "Positive" ~ "hilic"),
                  Metabolites = glue("{metabolite_name} ({column})")) %>%
    arrange(Metabolites) %>% 
    dplyr::mutate(Metabolites = factor(Metabolites, 
                               levels = .$Metabolites %>% unique()))  -> plotdata

  plotdata %>% 
    dplyr::select(any_of(var), Metabolites, value)  %>% 
    set_names(c("var", "Metabolites", "value")) %>% 
    dplyr::filter(!is.na(var)) -> a
  
    a %>%
      ggplot(aes(x = var, y = value)) +
      geom_point(color = "#7cadca", size = 2.5, alpha = 0.8) +                  # points
      geom_smooth(method = "lm", se = TRUE, color = "#b1564c",
                  fill = "#b1564c", alpha = 0.15, linewidth = 0.9) +            # fit + 95% CI
      stat_cor(method = "spearman",                                             # or "pearson"
             cor.coef.name = "rho",
             label.x.npc = "right", label.y.npc = "bottom", 
             hjust = 1,
             vjust = 0,
             size = 5) +
      theme( #plot parameters
            panel.border = element_rect(fill=NA, color="black", size=0.7, linetype="solid"),
            panel.background=element_blank(),
            axis.text = element_text(size = 14,  color = "black", face = "plain"),
            axis.title = element_text(size = 15, color = "black", face = "plain"), 
            strip.text = element_text(size = 14, color = "black", face = "plain")) +
        coord_cartesian() +
      scale_fill_manual(values=c("#7cadca", "#b1564c", "#E69F00")) + 
      scale_color_manual(values=c("#7cadca", "#b1564c", "#E69F00")) + 
      facet_wrap(~Metabolites , scales = "free_y", nrow = 5) +  
        labs(
          title = str_c("Distributions of metabolite levels across PMI "),
          y = "log2-transformed metabolite levels",
        x = var)
  
  ggsave(str_c(outpath, "/Visualizations/0_SM_figs", 
         "/Extended_Fig_metabolite level_comparison_across_[", var, "]_new", ".png"), 
         width = 16 ,
         height = 16,
         dpi = 600)

  print(str_c("Distributions of metabolite levels among group: ", var, " accomplished!"))

}

# execute 
purrr::map(c("Age_cat", "sex", "ApoE_cat", "ABC_cat", "Braak_cat", "CERAD_cat"),
     BoxPlot)
