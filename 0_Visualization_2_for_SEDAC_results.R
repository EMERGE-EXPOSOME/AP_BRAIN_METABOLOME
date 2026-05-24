rm(list=(ls()))
# 1- Prepare steps -------------------------------------
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
               "statmod", 
               "ggridges", 
               "kableExtra", 
               "RColorBrewer", # colors
               "ggtext", 
               "ggcorrplot", # correlation plot
               "patchwork", 
               "glue") # combine plots and control layouts

# set work directory
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# Create output folder
if(!file.exists(str_c(outpath,"/Visualizations"))){
  dir.create(str_c(outpath,"/Visualizations"))
}
if(!file.exists(str_c(outpath,"/Visualizations/SEDAC"))){
  dir.create(str_c(outpath,"/Visualizations/SEDAC"))
}

# 2-correlation plot -----------------------------
# load data
load(str_c(outpath, "/!Data_for_data_analysis_1.rdata"))

# df.c18 %>% 
#   dplyr::filter(!is.na(V1) & !is.na(SEDAC_PM25SO4_1y)) %>% 
#   .$SEDAC_O3_1y %>% summary()

df.c18 %>% 
  dplyr::filter(!is.na(V1) & !is.na(SEDAC_PM25SO4_1y)) %>% 
  dplyr::select(contains("SEDAC") & !ends_with("_15y")) -> tmp
colnames(tmp) %>% 
  as_tibble() %>% 
  dplyr::mutate(
    pollu = str_extract(value, "(?<=SEDAC_)[0-9a-zA-Z._]+(?=_1y|_3y|_5y|_10y)"), 
    window = str_extract(value, "(1|3|5|10)y$")
  )  -> tmp2
# extract unique factoral levels
unique(tmp2$pollu) -> pollu.lvl
c("1y", "3y", "5y", "10y") -> window.lvl
tmp2 %<>%
  dplyr::mutate(window = factor(window, levels = window.lvl), 
                pollu = factor(pollu, levels = pollu.lvl)) %>%
  arrange(pollu, window)
unique(tmp2$value) -> y.lvl
# re-arrange variables in tmp
tmp %<>% 
  dplyr::select(all_of(y.lvl))

# calculate correlation
cor.res = cor(tmp, use = "pairwise.complete.obs", method = "spearman")
cor.df = as.data.frame(as.table(cor.res))

####################Layout of correlation plot#########################
# column1: stacked blocks showing the time windows
# column2: stacked blocks showing the exposures
# column3: barplot(showing the VIF); stacked blocks showing the exposures and time windows; correlation matrix
#######################################################################

# column-3: barplot, stacked blocks, and correlation matrix ------------
# p3.4: correlation matrix
p3.4 = cor.df %>%
  filter(as.numeric(Var1) <= as.numeric(Var2)) %>%

  ggplot(aes(Var1, Var2, fill=Freq)) + 
  geom_tile(height=0.8, width=0.8) +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal() +
  scale_y_discrete(expand = expansion(c(0, 0))) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(
    axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(), 
    plot.margin = margin(0,0,0,0)
  )

# p3.3: stack blocks showing the exposures 
# assign factor levels
tmp2 %<>% 
  dplyr::mutate(
    value = factor(value, levels = y.lvl),
    pollu = factor(pollu, levels = pollu.lvl), 
    window = factor(window, levels = window.lvl))
# specify colors
color.pollu = c(
  "#BEB8DC", "#F28E2B", "#5D6572", "#8ECFC9", "#D37295",
  "#1F77B4", "#C4D6E5", "#9C755F", "#8A84BE", "#FF9DA7",
  "#4E79A7", "#EDC948", "#6C8EAD", "#B07AA1", "#76B7B2",
  "#A0C1D1", "#E15759", "#59A5D8")
color.window = c("#9dc0da", "#548dca", "#106caa", "#08179c")

p3.3 = tmp2 %>% 
  #dplyr::mutate(Pathways = factor(Pathways, levels = unique(.$Pathways))) %>% 
  ggplot(aes(x = value, y = 0)) + 
  geom_tile(aes(fill = pollu), alpha = 0.8) + 
  scale_fill_manual(values = color.pollu) + 
  scale_y_continuous(expand = expansion(c(0, 0))) +
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right", 
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank(),
        plot.margin   = margin(t = 0, r = 0, b = 0, l = 0)
      ) + 
  guides(color= "none",  size="none")

# p3.2: stack blocks showing the time windows
p3.2 = tmp2 %>% 
  #dplyr::mutate(Pathways = factor(Pathways, levels = unique(.$Pathways))) %>% 
  ggplot(aes(x = value, y = 0)) + 
  geom_tile(aes(fill = window), alpha = 0.8) + 
  scale_fill_manual(values = color.window) + 
  scale_y_continuous(expand = expansion(c(0, 0))) +
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right", 
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank(),
        plot.margin   = margin(t = 0, r = 0, b = 0, l = 0)
        ) + 
  guides(color= "none",  size="none") 

# p3.1: bar plot showing the VIF of each exposures before/after adding calendar year into models
vif = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables/VIF comparison between models.xlsx"), sheet = "cleaned")

p3.1 = vif %>% 
  dplyr::mutate(expo = factor(expo, levels = vif$expo)) %>% 
  pivot_longer(year_ex:year_in, names_to = "group", values_to = "value") %>% 
  dplyr::mutate(group = case_when(group=="year_ex" ~ "Exclude calendar year", 
                                  group=="year_in" ~ "Include calendar year", 
                                  TRUE ~ NA)) %>% 
  ggplot(aes(x = expo, y = log10(value), fill = group)) + 
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), expand = c(0, 0)) +
  scale_fill_manual(values = c("#f1e758", "#53a6ff")) +
  theme_bw() + 
  theme(
    # axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.x   = element_blank(),
    plot.margin   = margin(t = 1, r = 0, b = 0, l = 0),
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank())  + 
  labs(title = "Correlation and VIF of pollutants")

# merge column3
column3 = (p3.1/p3.2/p3.3/p3.4) + plot_layout(heights = c(0.4, 0.1, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right", 
        axis.text = element_text(size = 18))

# column2: stack plot for windows ----------------
# p2.4: stack blocks showing the time windows
p2.4 = tmp2 %>% 
  #dplyr::mutate(Pathways = factor(Pathways, levels = unique(.$Pathways))) %>% 
  ggplot(aes(x = 0, y = value)) + 
  geom_tile(aes(fill = window), alpha = 1) + 
  scale_fill_manual(values = color.window) + 
  scale_x_continuous(expand = expansion(c(0, 0))) +
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none", 
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank(),
        plot.margin   = margin(t = 0, r = 0, b = 0, l = 0)
        ) + 
  guides(color= "none",  size="none") 

# p2.1=p2.3: blank region
p2.1 = plot_spacer()
p2.2 = plot_spacer()
p2.3 = plot_spacer()

# merge column2
column2 = (p2.1/p2.2/p2.3/p2.4) + plot_layout(heights = c(0.4, 0.1, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "none", 
        axis.text = element_text(size = 18))

# column1: stack plot for pollutants ----------------
# p1.4: stack blocks showing the pollutants
p1.4 = tmp2 %>% 
  #dplyr::mutate(Pathways = factor(Pathways, levels = unique(.$Pathways))) %>% 
  ggplot(aes(x = 0, y = value)) + 
  geom_tile(aes(fill = pollu), alpha = 1) + 
  scale_fill_manual(values = color.pollu) + 
  scale_x_continuous(expand = expansion(c(0, 0))) +
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none", 
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank(),
        plot.margin   = margin(t = 0, r = 0, b = 0, l = 0)
        ) + 
  guides(color= "none",  size="none") 

# p1.1=p1.3: blank region
p1.1 = plot_spacer()
p1.2 = plot_spacer()
p1.3 = plot_spacer()

# merge column1
column1 = (p1.1/p1.2/p1.3/p1.4) + plot_layout(heights = c(0.4, 0.1, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "none", 
        axis.text = element_text(size = 18))

# merge column1+2+3 ------------------------------
(column1 | column2 | column3) + plot_layout(widths = c(0.1, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        # legend.position = "right", 
        axis.text = element_text(size = 15))

# save plot - with y axis tick values (log10 VIF) ----
ggsave(str_c(getwd(), "/0_Manuscript/figures/Correlation_VIF", ".png"), 
       height = 10, width = 10, dpi = 600)

# save plot without tick values -----------------------
p3.1 = vif %>% 
  dplyr::mutate(expo = factor(expo, levels = vif$expo)) %>% 
  pivot_longer(year_ex:year_in, names_to = "group", values_to = "value") %>% 
  dplyr::mutate(group = case_when(group=="year_ex" ~ "Exclude calendar year", 
                                  group=="year_in" ~ "Include calendar year", 
                                  TRUE ~ NA)) %>% 
  ggplot(aes(x = expo, y = log10(value), fill = group)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), expand = c(0, 0)) +
  theme_bw() + 
  theme(
    axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.x   = element_blank(),
    plot.margin   = margin(t = 1, r = 0, b = 0, l = 0),
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank())  + 
  labs(title = "Correlation and VIF of pollutants")

# merge
column3 = (p3.1/p3.2/p3.3/p3.4) + plot_layout(heights = c(0.4, 0.1, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right", 
        axis.text = element_text(size = 18))
(column1 | column2 | column3) + plot_layout(widths = c(0.1, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        # legend.position = "right", 
        axis.text = element_text(size = 15))

# save plot - with y axis tick values (log10 VIF) 
ggsave(str_c(getwd(), "//0_Manuscript/figures/Correlation_VIF_without_tick_value", ".png"), 
       height = 10, width = 10, dpi = 600)

######################################################################################
# UPDATE: only show 10-year in main text
cor.df %<>% 
  filter(str_detect(Var1, "_10y$") & str_detect(Var2, "_10y$"))

# column-3: barplot, stacked blocks, and correlation matrix ------------
# p3.3: correlation matrix
p3.3 = cor.df %>%
  filter(as.numeric(Var1) <= as.numeric(Var2)) %>%
  ggplot(aes(Var1, Var2, fill=Freq)) + 
  geom_tile(height=0.8, width=0.8) +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal() +
  scale_y_discrete(expand = expansion(c(0, 0))) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(
    axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(), 
    plot.margin = margin(0,0,0,0)
  )

# p3.2: stack blocks showing the exposures 
# assign factor levels
tmp2 %<>%
  filter(window == "10y") %>% 
  dplyr::mutate(
    value = factor(value, levels = y.lvl),
    pollu = factor(pollu, levels = pollu.lvl))
# specify colors
color.pollu = c(
  "#BEB8DC", "#F28E2B", "#5D6572", "#8ECFC9", "#D37295",
  "#1F77B4", "#C4D6E5", "#9C755F", "#8A84BE", "#FF9DA7",
  "#4E79A7", "#EDC948", "#6C8EAD", "#B07AA1", "#76B7B2",
  "#A0C1D1", "#E15759", "#59A5D8")

p3.2 = tmp2 %>% 
  #dplyr::mutate(Pathways = factor(Pathways, levels = unique(.$Pathways))) %>% 
  ggplot(aes(x = value, y = 0)) + 
  geom_tile(aes(fill = pollu), alpha = 0.8) + 
  scale_fill_manual(values = color.pollu) + 
  scale_y_continuous(expand = expansion(c(0, 0))) +
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right", 
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank(),
        plot.margin   = margin(t = 0, r = 0, b = 0, l = 0)
      ) + 
  guides(color= "none",  size="none")

# p3.1: bar plot showing the VIF of each exposures before/after adding calendar year into models
vif = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables/VIF comparison between models.xlsx"), sheet = "cleaned")
vif %<>% 
  filter(str_detect(expo, "_10y$"))
p3.1 = vif %>% 
  dplyr::mutate(expo = factor(expo, levels = vif$expo)) %>% 
  pivot_longer(year_ex:year_in, names_to = "group", values_to = "value") %>% 
  dplyr::mutate(group = case_when(group=="year_ex" ~ "Exclude calendar year", 
                                  group=="year_in" ~ "Include calendar year", 
                                  TRUE ~ NA)) %>% 
  ggplot(aes(x = expo, y = log10(value), fill = group)) + 
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), expand = c(0, 0)) +
  scale_fill_manual(values = c("#f1e758", "#53a6ff")) +
  theme_bw() + 
  theme(
    # axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.x   = element_blank(),
    plot.margin   = margin(t = 1, r = 0, b = 0, l = 0),
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank())  + 
  labs(title = "Correlation and VIF of pollutants")

# merge column2
column2 = (p3.1/p3.2/p3.3) + plot_layout(heights = c(0.4, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right", 
        axis.text = element_text(size = 18))

# column1: stack plot for pollutants ----------------
# p1.3: stack blocks showing the pollutants
p1.3 = tmp2 %>% 
  #dplyr::mutate(Pathways = factor(Pathways, levels = unique(.$Pathways))) %>% 
  ggplot(aes(x = 0, y = value)) + 
  geom_tile(aes(fill = pollu), alpha = 1) + 
  scale_fill_manual(values = color.pollu) + 
  scale_x_continuous(expand = expansion(c(0, 0))) +
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none", 
        axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.line.x   = element_blank(),
        plot.margin   = margin(t = 0, r = 0, b = 0, l = 0)
        ) + 
  guides(color= "none",  size="none") 

# p1.1=p1.2: blank region
p1.1 = plot_spacer()
p1.2 = plot_spacer()


# merge column1
column1 = (p1.1/p1.2/p1.3) + plot_layout(heights = c(0.4, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "none", 
        axis.text = element_text(size = 18))

# merge column1+2 ------------------------------
(column1 | column2) + plot_layout(widths = c(0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        # legend.position = "right", 
        axis.text = element_text(size = 15))

# save plot - with y axis tick values (log10 VIF) ----
ggsave(str_c(getwd(), "//0_Manuscript/figures/Correlation_VIF_10y_only", ".png"), 
       height = 10, width = 10, dpi = 600)

# save plot without tick values -----------------------
p3.1 = vif %>% 
  dplyr::mutate(expo = factor(expo, levels = vif$expo)) %>% 
  pivot_longer(year_ex:year_in, names_to = "group", values_to = "value") %>% 
  dplyr::mutate(group = case_when(group=="year_ex" ~ "Exclude calendar year", 
                                  group=="year_in" ~ "Include calendar year", 
                                  TRUE ~ NA)) %>% 
  ggplot(aes(x = expo, y = log10(value), fill = group)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), expand = c(0, 0)) +
  theme_bw() + 
  theme(
    axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.x   = element_blank(),
    plot.margin   = margin(t = 1, r = 0, b = 0, l = 0),
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank())  + 
  labs(title = "Correlation and VIF of pollutants")

# merge
column2 = (p3.1/p3.2/p3.3) + plot_layout(heights = c(0.4, 0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right", 
        axis.text = element_text(size = 18))
(column1 | column2) + plot_layout(widths = c(0.1, 5), guides = "collect") & 
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        # legend.position = "right", 
        axis.text = element_text(size = 15))

# save plot - with y axis tick values (log10 VIF) 
ggsave(str_c(getwd(), "//0_Manuscript/figures/Correlation_VIF_no_label_10y_only", ".png"), 
       height = 10, width = 10, dpi = 600)

#####################




# 3-Time trends of air pollutant exposure levels ---------
load(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
          "/!Data_for_DataAnalysis_2.rdata"))
# set work directory
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# see metabolic signatures that included in final MWAS analysis
df.c18 %>% 
  dplyr::select(matches("^V\\d+$")) %>% dim()
df.hilic %>% 
  dplyr::select(matches("^V\\d+$")) %>% dim()
# see how many participants included in MWAS analysis in different time windows
df.c18 %>% 
  dplyr::filter(!is.na(V15)) %>% nrow()
df.c18 %>% 
  dplyr::filter(!is.na(Rline_CO_5y) & !is.na(V1)) %>% nrow() # 102
df.c18 %>% 
  dplyr::filter(!is.na(SEDAC_PM25_10y) & !is.na(V1)) %>% nrow() # 114

# check out time trends of PM2.5
if(!file.exists(str_c(outpath,"/Visualizations/SEDAC/time_trends"))){
  dir.create(str_c(outpath,"/Visualizations/SEDAC/time_trends"))
}

# subset data
df.c18 %>% 
  dplyr::filter(!is.na(df.c18[["SEDAC_PM25NH4_1y"]]) & !is.na(V1)) -> tmp

# partcipant in brain tissue ADRD project
res = tibble(
  ID = tmp$GUID.x,
  death_yr = tmp$death_yr,
  Case_Number = tmp$Case_Number
)

# read annual air pollution data
fread(file = str_c(inpath, "/Air_pollution_data/Brain_exposure_matched_2025Feb18_from_Zhihao.csv")) -> ap.data
ap.data$year %>% range(na.rm = T) # 2000 2019

# for each participant
for(i in 1:nrow(res)){
  ap.data %>% 
    dplyr::filter(globalid == res$ID[i] & year <= res$death_yr[i]-1) %>% 
    dplyr::select(globalid, year, contains("SEDAC")) %>% 
    dplyr::rename(ID = "globalid") %>% 
    left_join(res[i, ], by = "ID") -> tmp2
    if(i ==1){
      plt_data = tmp2
    }else{
      plt_data = rbind(plt_data, tmp2)
    }
}

# dataset used for plot
plt_data %<>% 
  dplyr::select(year, contains("SEDAC"), death_yr, Case_Number) %>% 
  set_names(c("year", "NO2", "O3", "PM2.5", "NH4", "NO3", "OC", "EC", "SO4",
              "Bromine", "Calcium", "Copper", "Iron", "Lead", "Nickel", "Potassium", "Silicon", "Vanadium", "Zinc", 
              "death_yr", "subject"))
plt_data %>% 
  dplyr::select(NO2:Zinc) %>% colnames() -> pollu.lvl

# plotdata for time trend
plt_data %>% 
  pivot_longer(NO2:Zinc, names_to = "pollutant", values_to = "value") %>% 
  dplyr::mutate(pollutant = factor(pollutant, levels = pollu.lvl)) -> tmp
# extract mean value for each year
tmp %>% 
  group_by(pollutant, year) %>% 
  dplyr::summarise(mean = mean(value, na.rm = T)) %>% 
  dplyr::mutate(group = "mean")  ->  plt_data2 # data for mean line

ggplot(data = tmp, aes(x = year, y = value, group = subject)) + 
    geom_line(color = "#9b9898", alpha = 0.75) + 
    geom_line(data = plt_data2, aes(x = year, y = mean, group = group), 
              color = "#5351cf", alpha = 0.8, linewidth = 2) + 
    facet_wrap(pollutant~., scales = "free", ncol = 3) +
    labs(title = "Time trends of air pollutant exposure levels", 
         x = "Calendar year",
         y = str_c("Concentration")) +           
    theme_bw() + 
    theme(axis.text = element_text(size = 13), 
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold"))
ggsave(str_c(outpath,"/Visualizations/SEDAC/time_trends/", 
             "SM_fig.1_Time_trends",
                ".png"), width = 10, height = 16, dpi = 600)


# 4-Fig1: Ridgeline plot showing air pollutant distributions ---------
## 4.1-input ----
df.c18 %>% 
  dplyr::filter(!is.na(SEDAC_PM25NH4_1y) & !is.na(V1))  %>%
  dplyr::select(matches("SEDAC") & !matches("_15y$")) -> tmp_sedac
## 4.2-generate plotdata ----
tmp_sedac %>% 
  pivot_longer(everything(), names_to = "group", 
               values_to = "value") %>% 
  dplyr::mutate(windows = str_extract(group, "\\d+y$"), 
                pollutants = str_extract(group, "(?<=SEDAC_)\\w+(?=_\\d+y$)")) %>%
  dplyr::mutate(pollutants = str_replace_all(pollutants, "PM25", "PM2.5")) %>%
  dplyr::filter(str_detect(group, "_NO2_|_O3_|_PM25_")) %>% 
rbind(
tmp_sedac %>% 
  pivot_longer(everything(), names_to = "group", 
               values_to = "value") %>% 
  dplyr::mutate(windows = str_extract(group, "\\d+y$"), 
                pollutants = str_extract(group, "(?<=SEDAC_PM25)\\w+(?=_\\d+y$)")) %>% 
  dplyr::filter(str_detect(group, "_PM25\\w+_"))
) %>% 
rbind(
tmp_sedac %>% 
  pivot_longer(everything(), names_to = "group", 
               values_to = "value") %>% 
  dplyr::mutate(windows = str_extract(group, "\\d+y$"), 
                pollutants = str_extract(group, "(?<=SEDAC_PM2.5_)\\w+(?=_\\d+y$)")) %>% 
  dplyr::filter(str_detect(group, "_PM2.5_"))
) %>% 
  dplyr::select(-group) -> tmp
tmp$pollutants %>% unique()
# assign levels
tmp %<>% 
  dplyr::mutate(windows = factor(windows, levels = c("1y", "3y", "5y", "10y")), 
                pollutants = factor(pollutants, levels = pollu.lvl))

# plot
tmp %>%
    ggplot(aes(x = log2(value), y = windows, fill = windows)) +
      geom_density_ridges(alpha = 0.2, scale = 2, stat = "binline", bins = 20, color = "#9b9b9e", fill = "#9b9b9e") +
      geom_density_ridges(alpha = 0.6, scale = 2) +
      theme_bw() + 
      facet_wrap(pollutants~., scales = "free", ncol = 3) +
      scale_y_discrete(expand = c(0.01, 0.01)) + 
      scale_x_continuous(expand = c(0.01, 0.01)) +
      scale_fill_cyclical(values = c("#bbbcee", "#7f81f0", "#5659f0", "#090df0")) +
      labs(title = str_c("Distribution of log2-transformed air pollutant concentrations across time periods"),
           x = "log2(concentration)",
           y = "Frequencies") +
      theme(legend.position = "none", 
            axis.text.y = element_text(size = 16), 
            axis.text.x = element_text(size = 12), 
            axis.title = element_text(size = 20),
            strip.text = element_text(size = 16),
            plot.title = element_text(size = 20, face = "bold"))
ggsave(str_c(outpath, "/Visualizations/SEDAC/SM_fig.2_exposure_distri.png"), 
         height = 16, width = 10, dpi = 600)


# 3-Fig2: Hotplot showing significant metabolic signatures ---------
## 3.1-read data ----
### 3.1.1-main components and mixture fdr0.2----
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05", 
      "/Pvalue_counts_XtoM_SEDAC.csv")) %>% 
  group_by(Exposure) %>% 
  dplyr::filter(str_detect(Exposure, "_PM25_|_NO2_|_O3_")) %>% # include main components
  dplyr::mutate(rowSum = sum(P_fdr_0.2)) %>% 
  dplyr::select(Exposure, Type, rowSum) %>% 
  unique() -> fig2.p1
fig2.p1 %<>% 
  dplyr::mutate(pollutants = str_extract(Exposure, "(?<=SEDAC_)\\w+(?=_)"), 
                windows = str_extract(Exposure, "(?<=_)\\d+y$"))
# mixture results
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture/", 
      "/Pvalue_counts_XtoM_SEDAC.main_mixture.csv")) %>% 
  group_by(Exposure) %>% 
  dplyr::mutate(rowSum = sum(P_fdr_0.2)) %>% 
  dplyr::select(Exposure, Type, rowSum)  %>% 
  unique() -> fig2.p2
fig2.p2 %<>% 
  dplyr::mutate(pollutants = "mixture", 
                windows = str_extract(Exposure, "(?<=mix_)\\d+y$"))
fig2.fdr0.2 <- rbind(fig2.p1, fig2.p2) %>% 
  dplyr::mutate(group = "fdr_0.2")

### 3.1.1-main components and mixture raw0.05----
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05", 
      "/Pvalue_counts_XtoM_SEDAC.csv")) %>% 
  group_by(Exposure) %>% 
  dplyr::filter(str_detect(Exposure, "_PM25_|_NO2_|_O3_")) %>% # include main components
  dplyr::mutate(rowSum = sum(P_raw_0.05)) %>% 
  dplyr::select(Exposure, Type, rowSum) %>% 
  unique() -> fig2.p1
fig2.p1 %<>% 
  dplyr::mutate(pollutants = str_extract(Exposure, "(?<=SEDAC_)\\w+(?=_)"), 
                windows = str_extract(Exposure, "(?<=_)\\d+y$"))
# mixture results
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture", 
      "/Pvalue_counts_XtoM_SEDAC.main_mixture.csv")) %>% 
  group_by(Exposure) %>% 
  dplyr::mutate(rowSum = sum(P_raw_0.05)) %>% 
  dplyr::select(Exposure, Type, rowSum)  %>% 
  unique() -> fig2.p2
fig2.p2 %<>% 
  dplyr::mutate(pollutants = "mixture", 
                windows = str_extract(Exposure, "(?<=mix_)\\d+y$"))
fig2.raw0.05 <- rbind(fig2.p1, fig2.p2) %>% 
  dplyr::mutate(group = "raw_0.05")
# SEDAC main components
fig2.main.fix = rbind(fig2.fdr0.2, fig2.raw0.05)
# View(fig2.main.fix)
fig2.main.fix %>% 
  dplyr::mutate(pollutants = factor(pollutants, levels = .$pollutants %>% unique()), 
                windows = factor(windows, levels = c("1y", "3y", "5y", "10y"))) %>% 
  ggplot(aes(pollutants, windows)) +
  geom_tile(aes(fill = rowSum)) + 
  scale_fill_gradientn(colors = c("#bebebe7c", "#c0c0f5ff", "blue"), 
                       values = scales::rescale(c(0, 1, max(fig2.main.fix$rowSum)))) + 
  facet_grid(.~group, scales = "free", space = "free") + 
  theme_bw() + 
  theme(panel.spacing = unit(0, "cm", data = NULL), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        plot.title = element_text(size = 14, face = "bold")) + 
  labs(title = "Numbers of significant metabolic signatures")


### 3.1.2-PM2.5 + metal components and mixtures fdr0.2----
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05", 
      "/Pvalue_counts_XtoM_SEDAC.csv")) %>% 
  group_by(Exposure) %>% 
  dplyr::filter(!str_detect(Exposure, "_PM25_|_NO2_|_O3_")) %>% # include main components
  dplyr::mutate(rowSum = sum(P_fdr_0.2)) %>% 
  dplyr::select(Exposure, Type, rowSum) %>% 
  unique() -> fig2.p1
fig2.p1 %<>% 
  dplyr::mutate(pollutants = str_extract(Exposure, "(?<=(SEDAC_PM25|SEDAC_PM2.5_))\\w+(?=_\\d+y)"), 
                windows = str_extract(Exposure, "(?<=_)\\d+y$"))
# mixture results
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture", 
      "/Pvalue_counts_XtoM_SEDAC.pm25_mixture.csv")) %>% 
rbind(
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture", 
      "/Pvalue_counts_XtoM_SEDAC.metal_mixture.csv"))  
) %>% 
rbind(
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture", 
      "/Pvalue_counts_XtoM_SEDAC.pm25.full_mixture.csv"))  
) %>% 
  group_by(Exposure, Type) %>% 
  dplyr::mutate(rowSum = sum(P_fdr_0.2)) %>% 
  dplyr::select(Exposure, Type, rowSum)  %>% 
  unique() -> fig2.p2
fig2.p2 %<>% 
  dplyr::mutate(pollutants = str_c("mixture_", str_extract(Type, "(?<=SEDAC.)[a-zA-Z0-9.]+")), 
                windows = str_extract(Exposure, "(?<=mix_)\\d+y$"))
fig2.fdr0.2 <- rbind(fig2.p1, fig2.p2) %>% 
  dplyr::mutate(group = "fdr_0.2")

### 3.1.3-PM2.5 + metal components and mixtures raw0.05----
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05", 
      "/Pvalue_counts_XtoM_SEDAC.csv")) %>% 
  group_by(Exposure) %>% 
  dplyr::filter(!str_detect(Exposure, "_PM25_|_NO2_|_O3_")) %>% # include main components
  dplyr::mutate(rowSum = sum(P_raw_0.05)) %>% 
  dplyr::select(Exposure, Type, rowSum) %>% 
  unique() -> fig2.p1
fig2.p1 %<>% 
  dplyr::mutate(pollutants = str_extract(Exposure, "(?<=(SEDAC_PM25|SEDAC_PM2.5_))\\w+(?=_\\d+y)"), 
                windows = str_extract(Exposure, "(?<=_)\\d+y$"))
# mixture results
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture", 
      "/Pvalue_counts_XtoM_SEDAC.pm25_mixture.csv")) %>% 
rbind(
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture", 
      "/Pvalue_counts_XtoM_SEDAC.metal_mixture.csv"))  
) %>% 
rbind(
fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture", 
      "/Pvalue_counts_XtoM_SEDAC.pm25.full_mixture.csv"))  
) %>% 
  group_by(Exposure, Type) %>% 
  dplyr::mutate(rowSum = sum(P_raw_0.05)) %>% 
  dplyr::select(Exposure, Type, rowSum)  %>% 
  unique() -> fig2.p2
fig2.p2 %<>% 
  dplyr::mutate(pollutants = str_c("mixture_", str_extract(Type, "(?<=SEDAC.)[a-zA-Z0-9.]+")), 
                windows = str_extract(Exposure, "(?<=mix_)\\d+y$"))
fig2.raw0.05 <- rbind(fig2.p1, fig2.p2) %>% 
  dplyr::mutate(group = "raw_0.05")
fig2.pm25 = rbind(fig2.fdr0.2, fig2.raw0.05)

# merge 
fig2 = rbind(fig2.main.fix, fig2.pm25) %>% 
  dplyr::mutate(pollutants = case_when(pollutants == "PM25" ~ "PM2.5",
                                       pollutants == "mixture" ~ "Mixture(NO2+O3+PM2.5)",
                                       pollutants == "mixture_pm25" ~ "Mixture(PM2.5 species)",
                                       pollutants == "mixture_metal" ~ "Mixture(metal components)",
                                       pollutants == "mixture_pm25.full" ~ "Mixture(PM2.5 species + metals)",
                                       TRUE ~ pollutants), 
                group = case_when(group == "fdr_0.2" ~ "FDR q<0.2",
                                  group == "raw_0.05" ~ "Raw p<0.05", 
                                  TRUE ~ group))

# assign levels
fig2$pollutants %>% unique()  -> pollu.lvl
new.lvl = pollu.lvl[c(1:9,20,10:19,21,22)]

## 3.2-draw plots --------------------
fig2 %>%
  dplyr::mutate(pollutants = factor(pollutants, levels = rev(new.lvl)), 
                windows = factor(windows, levels = c("1y", "3y", "5y", "10y"))) %>% 
  ggplot(aes(windows, pollutants)) +
  geom_tile(aes(fill = rowSum), alpha=1) + 
  scale_fill_gradientn(colors = c("#ddd9d97c", "#eaedf8", "blue"), 
                       values = scales::rescale(c(0, 1, max(fig2$rowSum))), 
                       limits = c(0, max(fig2$rowSum))) + 
  facet_grid(.~group, scales = "free", space = "free") + 
  theme_bw() + 
  theme(panel.spacing = unit(0, "cm", data = NULL), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold")) + 

  labs(title = "Numbers of significant metabolic signatures")
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig2_MWAS", ".png"), 
       height = 12, width = 10, dpi = 600)
# export table
fig2 %>% 
  ungroup() %>% 
  dplyr::select(-Exposure) %>% 
  pivot_wider(names_from = "windows", values_from = "rowSum") %>% 
  fwrite(str_c(outpath, "/Visualizations/SEDAC/Fig2_MWAS_summary.csv"))


df = readxl::read_xlsx(str_c(work.dir, "/0_Manuscript/tables/2_Sample_size_comparison.xlsx"))


# merge 
df %<>% 
  pivot_longer('1y':'10y', names_to = "windows", values_to = "values")

# assign levels
df$Pollutant %>% unique()  -> pollu.lvl
new.lvl = pollu.lvl[c(1:9,20,10:19,21,22)]

## 3.2-draw plots --------------------
df %>%
  dplyr::mutate(pollutants = factor(Pollutant, levels = rev(new.lvl)), 
                windows = factor(windows, levels = c("1y", "3y", "5y", "10y"))) %>% 
  ggplot(aes(windows, pollutants)) +
  geom_tile(aes(fill = values), alpha=1) + 
  scale_fill_gradientn(colors = c("#ddd9d97c", "#eaedf8", "blue"), 
                       values = scales::rescale(c(0, 1, max(df$values))), 
                       limits = c(0, max(df$values))) + 
  facet_grid(.~Group, scales = "free", space = "free") + 
  theme_bw() + 
  theme(panel.spacing = unit(0, "cm", data = NULL), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold")) + 

  labs(title = "Numbers of significant metabolic signatures")
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig2_MWAS", ".png"), 
       height = 12, width = 10, dpi = 600)




# 4-UpSetR showing intersection sets of metabolic signature numbers ---------
library(visreg)
library(UpSetR)
library(ComplexUpset)
library(ggbreak)

## 4.1-read data ----
sedac = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
        "/MWAS_TRAP_Metabolome_SEDAC.rds"))
mix = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
        "/MWAS_TRAP_Metabolome_mixture.rds"))
mix$sedac.main$c18$coef[["5y"]]


# function to extract data
F_UpSet.Data = function(type="SEDAC"){
  switch(type, 
    "RLINE"={
      lst = rline # individual pollutants
      lst.2 = mix$rline

    }, 
    "SEDAC"={
      lst = sedac # individual pollutants
      lst.2.main = mix$sedac.main
      lst.2.pm25 = mix$sedac.pm25
      lst.2.metal = mix$sedac.metal
      lst.2.full = mix$sedac.pm25.full
    })
  
  list.names(lst) -> list.names
  purrr::map(list.names, 
    ~{
      lst[[.x]]$c18 %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
      dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when raw p<0.001, marked as 1, otherwise as 0
                    Components = str_replace(Components, "V", "c")) %>% 
      dplyr::select(Components, Pvalue_a) %>%
      rbind(
      lst[[.x]]$hilic %>% 
        dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH")) %>% 
        dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when raw p<0.001, marked as 1, otherwise as 0
                      Components = str_replace(Components, "V", "h")) %>% 
        dplyr::select(Components, Pvalue_a)
      ) %>%
      dplyr::mutate(group = .x, 
                    type  = type) %>% 
      dplyr::select(type, Components, group, everything())
    })  %>% 
    rbindlist() %>% 
    dplyr::mutate(group = str_replace_all(group, "SEDAC_|PM2.5_", "")) %>% 
    dplyr::mutate(group = case_when(str_detect(group, "PM25[a-zA-Z0-9]+_\\d+y$") ~ str_extract(group, "(?<=PM25)[a-zA-Z0-9]+_\\d+y$"), 
                                    TRUE ~ group)) %>% 
    pivot_wider(names_from = "group", values_from = "Pvalue_a") -> tmp

  if(type == "SEDAC"){
    list.names.2 = c("1y", "3y", "5y", "10y")
    # main pollutants
    purrr::map(list.names.2, 
      ~{
        lst.2.main$c18$coef[[.x]] %>% 
        dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
        dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                      Components = str_replace(Outcome, "V", "c")) %>% 
        dplyr::select(Components, Pvalue_a) %>%
        rbind(
        lst.2.main$hilic$coef[[.x]] %>% 
          dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
          dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                        Components = str_replace(Outcome, "V", "h")) %>% 
          dplyr::select(Components, Pvalue_a)
        ) %>%
        dplyr::mutate(group = str_c("Mixture(NO2+O3+PM2.5)_", .x), 
                      type  = type) %>% 
        dplyr::select(type, Components, group, everything())
      })  %>% 
      rbindlist()  -> tmp.p1
    
    # pm2.5 species
    purrr::map(list.names.2, 
      ~{
        lst.2.pm25$c18$coef[[.x]] %>% 
        dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
        dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                      Components = str_replace(Outcome, "V", "c")) %>% 
        dplyr::select(Components, Pvalue_a) %>%
        rbind(
        lst.2.pm25$hilic$coef[[.x]] %>% 
          dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
          dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                        Components = str_replace(Outcome, "V", "h")) %>% 
          dplyr::select(Components, Pvalue_a)
        ) %>%
        dplyr::mutate(group = str_c("Mixture(PM2.5 species)_", .x), 
                      type  = type) %>% 
        dplyr::select(type, Components, group, everything())
      })  %>% 
      rbindlist()  -> tmp.p2

    # pm2.5 metal components
    purrr::map(list.names.2, 
      ~{
        lst.2.metal$c18$coef[[.x]] %>% 
        dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
        dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                      Components = str_replace(Outcome, "V", "c")) %>% 
        dplyr::select(Components, Pvalue_a) %>%
        rbind(
        lst.2.metal$hilic$coef[[.x]] %>% 
          dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
          dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                        Components = str_replace(Outcome, "V", "h")) %>% 
          dplyr::select(Components, Pvalue_a)
        ) %>%
        dplyr::mutate(group = str_c("Mixture(metal components)_", .x), 
                      type  = type) %>% 
        dplyr::select(type, Components, group, everything())
      })  %>% 
      rbindlist()  -> tmp.p3

    # pm2.5 full components
    purrr::map(list.names.2, 
      ~{
        lst.2.full$c18$coef[[.x]] %>% 
        dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
        dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                      Components = str_replace(Outcome, "V", "c")) %>% 
        dplyr::select(Components, Pvalue_a) %>%
        rbind(
        lst.2.full$hilic$coef[[.x]] %>% 
          dplyr::mutate(P_fdr = p.adjust(P.value, method="BH")) %>% 
          dplyr::mutate(Pvalue_a = if_else(P_fdr<0.2, 1, 0), # when FDR<0.2, marked as 1, otherwise as 0
                        Components = str_replace(Outcome, "V", "h")) %>% 
          dplyr::select(Components, Pvalue_a)
        ) %>%
        dplyr::mutate(group = str_c("Mixture(PM2.5 species + metals)_", .x), 
                      type  = type) %>% 
        dplyr::select(type, Components, group, everything())
      })  %>% 
      rbindlist()  -> tmp.p4
    
    tmp.2 = rbind(tmp.p1, tmp.p2, tmp.p3, tmp.p4) %>% 
      pivot_wider(names_from = "group", values_from = "Pvalue_a")
    
    tmp %>% 
      full_join(tmp.2, by = c("type", "Components")) -> tmp
  }# end if type==sedac
  
  print(str_c("UpSet data for visit ", type, " has been exrtacted !"))
}

# rbind
tmp %>% 
  as_tibble() -> df

## 4.3 draw Fig -----------
### 4.3.1  for SEDAC-----
F_UpSet.Plot = function(source = "pollutant",
                        window = "1y"){
  switch(source,
         "pollutant" = df <- tmp %>% dplyr::select(type, Components, contains(c("PM25_", "NO2_", "O3_", "Mixture(NO2+O3+PM2.5)"))),
         "pm.species" = df <- tmp %>% dplyr::select(type, Components, contains(c("NH4_", "NO3_", "OC_", "EC_", "SO4_", "Mixture(PM2.5 species)"))),
         "pm.components" = df <- tmp %>% dplyr::select(type, Components, contains(c("Bromine_","Calcium_","Copper_","Iron_","Lead_","Nickel_","Potassium_","Silicon_","Vanadium_","Zinc_","Mixture(metal components)_"))),
         "manual" = df <- tmp %>% dplyr::select(type, Components, contains(c("Vanadium_", "Silicon_", "SO4_", "NH4_", "O3_", "Mixture(NO2+O3+PM2.5)")))
        )

  df %>% 
    dplyr::select(type, Components, contains(window)) %>% 
    dplyr::mutate(column = if_else(str_detect(Components, "^c"), "c18", "hilic")) -> fig

  fig %>% dplyr::select(contains(window)) %>% colnames() -> grp

  # UpSet plot
  upset(fig, 
        grp, 
        name='TRAP components', 
        min_degree= 1,
        height_ratio = 0.2,
        width_ratio=0.1,
        base_annotations=list(
          'Intersection metabolite numbers across TRAP components'=intersection_size(
              counts=FALSE,
              mapping=aes(fill=column)
          )
          + scale_fill_manual(values=c(
              'c18'='#8b4446', 'hilic'='#487da8'
          ))
          + scale_y_continuous(expand = c(0, 0))
          + labs(title = str_c(source," data with ", window, " time window"))
          #+ scale_y_break(c(100, 1900))
        ),
        matrix = (
          intersection_matrix(
              segment=geom_segment(
                  # linetype='dotted'
                  color = "#302e2e"
              ),
              outline_color=list(
                  active='#302e2e',
                  inactive='#929090'
              ),
              geom=geom_point(
                  # shape='square',
                  size=3
              )
          ) # intersection_matrix end
          + scale_color_manual(
                values=c('TRUE'='#302e2e', 'FALSE'='#929090'),
                labels=c('TRUE'='yes', 'FALSE'='no')
            )
          + scale_y_discrete(expand = c(0.1, 0.1))
        ), # matrix end
        stripes='white',
        set_sizes=FALSE) 
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig3_",source, "_", window,".png"), 
         height = 11, width = 7, dpi = 600)
}

# main pollutant
purrr::map2(rep("pollutant", 4), 
            c("1y", "3y", "5y", "10y"), 
            F_UpSet.Plot)
# pm species, such as OC, EC
purrr::map2(rep("pm.species", 4), 
            c("1y", "3y", "5y", "10y"), 
            F_UpSet.Plot)
# pm components, such as metal elements
purrr::map2(rep("pm.components", 4), 
            c("1y", "3y", "5y", "10y"), 
            F_UpSet.Plot)
# manual
purrr::map2(rep("manual", 4), 
            c("1y", "3y", "5y", "10y"), 
            F_UpSet.Plot)


## 5 draw venn plots ----
### 5.1 betwen visits ----
# read data-individual SEDAC POLLUTANT AND PM SPECIES & COMPONENTS
res1 = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
                       "/MWAS_TRAP_Metabolome_SEDAC.rds"))
res1 %>% list.names()
# read data-mixture results
res2 = readRDS(str_c(outpath, "/DataAnalysis_1_TRAP_Metabolome/Sensitive_analyses_using maximum datasets", 
                       "/MWAS_TRAP_Metabolome_mixture.rds"))
res2 %>% list.names()
res1$SEDAC_NO2_1y
res2$sedac.pm25.full$c18$coef %>% list.names()
# res.mix$v1$c18$coef %>% list.names()
# res.mix$v1$hilic$coef %>% list.names()
# res.mix$v2$c18$coef %>% list.names()
# res.mix$v2$hilic$coef %>% list.names()

# integrate data into one single dataframe
#c18
map2(res1, res1 %>% list.names(), 
     ~{.x$c18 %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_extract(.y, "(?<=SEDAC_)[a-zA-Z0-9._]+"), 
                    Pollutant = str_extract(Exposure, "[a-zA-Z0-9._]+(?=_(1|3|5|10)y$)"), 
                    Window = str_extract(Exposure, "(1|3|5|10)y$")) %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> c18.p1 # individual pollutants
map2(res2$sedac.main$c18$coef, res2$sedac.main$c18$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("TRAP mix_", .y), 
                    Pollutant = "TRAP mix") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> c18.p2 # mixture exposure-TRAP components (PM2.5, O3, NO2)
map2(res2$sedac.pm25$c18$coef, res2$sedac.pm25$c18$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("PM2.5 species_", .y), 
                    Pollutant = "PM2.5 species") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> c18.p3 # mixture exposure-PM2.5 species
map2(res2$sedac.metal$c18$coef, res2$sedac.metal$c18$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("PM2.5 metal_", .y), 
                    Pollutant = "PM2.5 metal") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> c18.p4 # mixture exposure-pm25 metal components
map2(res2$sedac.pm25.full$c18$coef, res2$sedac.pm25.full$c18$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("PM2.5 species and metals_", .y), 
                    Pollutant = "PM2.5 species and metals") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> c18.p5 # mixture exposure-pm25 species plus components
# merge
c18 = rbind(c18.p1, c18.p2, c18.p3, c18.p4, c18.p5)

#hilic
map2(res1, res1 %>% list.names(), 
     ~{.x$hilic %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_extract(.y, "(?<=SEDAC_)[a-zA-Z0-9._]+"), 
                    Pollutant = str_extract(Exposure, "[a-zA-Z0-9._]+(?=_(1|3|5|10)y$)"), 
                    Window = str_extract(Exposure, "(1|3|5|10)y$")) %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> hilic.p1 # individual pollutants
map2(res2$sedac.main$hilic$coef, res2$sedac.main$hilic$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("TRAP mix_", .y), 
                    Pollutant = "TRAP mix") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> hilic.p2 # mixture exposure-TRAP components (PM2.5, O3, NO2)
map2(res2$sedac.pm25$hilic$coef, res2$sedac.pm25$hilic$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("PM2.5 species_", .y), 
                    Pollutant = "PM2.5 species") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> hilic.p3 # mixture exposure-PM2.5 species
map2(res2$sedac.metal$hilic$coef, res2$sedac.metal$hilic$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("PM2.5 metal_", .y), 
                    Pollutant = "PM2.5 metal") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> hilic.p4 # mixture exposure-pm25 metal components
map2(res2$sedac.pm25.full$hilic$coef, res2$sedac.pm25.full$hilic$coef %>% list.names(), 
     ~{.x %>% 
      dplyr::rename(Pvalue_a = "P.value", 
                    Window = "Windows", 
                    Components = "Outcome") %>% 
      dplyr::mutate(P_fdr = p.adjust(Pvalue_a, method="BH"), 
                    Exposure = str_c("PM2.5 species and metals_", .y), 
                    Pollutant = "PM2.5 species and metals") %>% 
      dplyr::select(Exposure, Pollutant, Window, P_fdr, Pvalue_a, Components)
      }
    ) %>% 
rbindlist()  -> hilic.p5 # mixture exposure-pm25 species plus components
# merge
hilic = rbind(hilic.p1, hilic.p2, hilic.p3, hilic.p4, hilic.p5)

# modify pollutant names
c18$Pollutant %>% unique() -> pollu.lst
pollu.lst
# 
c18 %<>% 
  dplyr::mutate(Pollutant = case_when(Pollutant %in% pollu.lst[1:3] ~ str_c("1_", Pollutant), 
                                      Pollutant %in% pollu.lst[4:8] ~ str_c("2_", str_replace_all(Pollutant, "PM25", "")), 
                                      Pollutant %in% pollu.lst[9:18] ~ str_c("3_", str_replace_all(Pollutant, "PM2.5_", "")), 
                                      Pollutant == "TRAP mix" ~ "4_TRAP mix", 
                                      Pollutant == "PM2.5 species" ~ "5_PM2.5 species mix", 
                                      Pollutant == "PM2.5 metal" ~ "6_PM2.5 metal mix", 
                                      Pollutant == "PM2.5 species and metals" ~ "7_PM2.5 species and metals mix",
                                      TRUE ~ Pollutant))
# 
hilic %<>% 
  dplyr::mutate(Pollutant = case_when(Pollutant %in% pollu.lst[1:3] ~ str_c("1_", Pollutant), 
                                      Pollutant %in% pollu.lst[4:8] ~ str_c("2_", str_replace_all(Pollutant, "PM25", "")), 
                                      Pollutant %in% pollu.lst[9:18] ~ str_c("3_", str_replace_all(Pollutant, "PM2.5_", "")), 
                                      Pollutant == "TRAP mix" ~ "4_TRAP mix", 
                                      Pollutant == "PM2.5 species" ~ "5_PM2.5 species mix", 
                                      Pollutant == "PM2.5 metal" ~ "6_PM2.5 metal mix", 
                                      Pollutant == "PM2.5 species and metals" ~ "7_PM2.5 species and metals mix",
                                      TRUE ~ Pollutant))

### 5.2 betwen different windows, for PM2.5 specifically ----
# devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)
# we could also add mixture exposure results in SM file...but basically, results would be the similar as those in PM2.5 settings.
Venn_TW = function(pollu = "PM25", 
                   threshold = "Raw p<0.05") {

  switch(threshold,
    "FDR q<0.2" = {
        c18 %>% 
          dplyr::filter(str_detect(Pollutant, pollu)) %>% 
          dplyr::filter(P_fdr < 0.2) %>% 
          dplyr::mutate(Components = str_replace(Components, "V", "c")) %>% 
        rbind(
        hilic %>% 
          dplyr::filter(str_detect(Pollutant, pollu)) %>% 
          dplyr::filter(P_fdr < 0.2) %>% 
          dplyr::mutate(Components = str_replace(Components, "V", "h"))
        )  -> venn

        color = "#eb2d2df3"
        pvalue = "fdr0.2"
      },
    "Raw p<0.05" = {
        c18 %>% 
          dplyr::filter(str_detect(Pollutant, pollu)) %>% 
          dplyr::filter(Pvalue_a < 0.05) %>% 
          dplyr::mutate(Components = str_replace(Components, "V", "c")) %>% 
        rbind(
        hilic %>% 
          dplyr::filter(str_detect(Pollutant, pollu)) %>% 
          dplyr::filter(Pvalue_a < 0.05) %>% 
          dplyr::mutate(Components = str_replace(Components, "V", "h"))
        )  -> venn

        color = "#5454ebed"
        pvalue = "raw0.05"
      }
  )
  
  #construct list for venn plot
  venn.list = list(
    "1y" = venn %>% dplyr::filter(Pollutant == pollu & Window == "1y") %>% 
      .$Components %>% unique(),
    "3y" = venn %>% dplyr::filter(Pollutant == pollu & Window == "3y") %>% 
      .$Components %>% unique(),
    "5y" = venn %>% dplyr::filter(Pollutant == pollu & Window == "5y") %>% 
      .$Components %>% unique(),
    "10y" = venn %>% dplyr::filter(Pollutant == pollu & Window == "10y") %>% 
      .$Components %>% unique()
  )
  names(venn.list) = str_c(pollu, "_", c("1y", "3y", "5y", "10y"))

  # plot
  ggVennDiagram(venn.list, label_alpha=0, label = "count",label_size = 8) + 
    scale_fill_gradient(low="#f0dddd",high =color) + 
    theme()
  ggsave(str_c(outpath, "/Visualizations/SEDAC/", pollu,  "_", pvalue, ".png"), 
         height = 10, width = 10.5, dpi = 600)
}

# list of pollutants: 
# raw
purrr::pmap(list(c18$Pollutant %>% unique(), 
                 rep("Raw p<0.05", times = length(c18$Pollutant %>% unique()))), 
            Venn_TW)
# fdr
purrr::pmap(list(c18$Pollutant %>% unique(), 
                 rep("FDR q<0.2", times = length(c18$Pollutant %>% unique()))), 
            Venn_TW)


####################################################################################


## 6-1:Pathway analysis viz- SEDAC (raw signal p<0.05 as input)-----------------------------------
file_names = list.files(path = str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05"), pattern= "metapone_pathways_SEDAC_[a-zA-Z0-9._]+.xlsx$")
expo_names = file_names %>% str_extract("(?<=metapone_pathways_SEDAC_)[a-zA-Z0-9._]+(?=.xlsx)")
map2(file_names, expo_names, 
  ~readxl::read_xlsx(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/", .x)) %>% 
    dplyr::mutate(Exposure = .y) %>% 
    dplyr::mutate(Expo_grp = str_extract(Exposure, "[a-zA-Z0-9._]+(?=_\\d+y$)"),
                  Windows = str_extract(Exposure, "(?<=_)\\d+y$")) %>% 
    dplyr::select(Exposure, Expo_grp, Windows, everything())) %>% 
  rbindlist() -> plotdata

file_names = list.files(path = str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture"), pattern= "metapone_pathways_mixture_sedac.[a-zA-Z0-9._]+.xlsx$")
expo_names = file_names %>% str_extract("(?<=metapone_pathways_mixture_sedac.)[a-zA-Z0-9._]+(?=.xlsx)")
map2(file_names, expo_names, 
  ~readxl::read_xlsx(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05/TRAP_mixture/", .x)) %>% 
    dplyr::mutate(Exposure = str_c("mixture_", .y)) %>% 
    dplyr::mutate(Expo_grp = str_extract(Exposure, "[a-zA-Z0-9._]+(?=_\\d+y$)"),
                  Windows = str_extract(.y, "(?<=_)\\w+$")) %>% 
    dplyr::select(Exposure, Expo_grp, Windows, everything())) %>% 
  rbindlist() -> plotdata2

plt.data = rbind(plotdata, plotdata2)
# plt.data %>% View()

# unique pathways, there were a total of 230 pathways identified among all group settings 
plt.data$Pathways %>% unique() %>% length()

###Plot----------------------------
#deal with p=0
plt.data$P <- ifelse(plt.data$p_value == 0, 0.0000001, plt.data$p_value)
plt.data %>% colnames()
# raw p<0.05
plt.data %>% 
  dplyr::mutate(P = if_else(p_value == 0, 0.0000001, p_value))  -> tmp

# export for manually checking pathway names
tmp %>% 
  fwrite(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05", 
               "/0_pathway_analysis_results.csv"))


#####################################################################
# read modified tmp back to environment
tmp = readxl::read_xlsx(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.05", 
               "/0_pathway_analysis_results_raw0.05[manually_modified].xlsx"))
#####################################################################
library(patchwork)

tmp %>% 
  dplyr::arrange(ptw.main.category) %>% 
  .$ptw.main.category %>% 
  unique() -> ptw.cate.lvl # pathway main category levels, n=23

tmp %>% 
  dplyr::select(Pathways, ptw.main.category) %>% 
  unique() %>% 
  arrange(ptw.main.category) %>% 
  .$Pathways  -> ptw.lvl # pathway levels

ptw.color23 = c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
  "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
  "#393B79", "#637939", "#8C6D31"
)
expo.color22 = c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#97a884",
  "#E6AB02", "#A6761D", "#a19292", "#A6CEE3", "#535551",
  "#FDBF6F", "#CAB2D6"
)

# rank pathways by: 1) pathway main category; and 2) alphabet name order
colnames(tmp)
tmp %<>% 
  dplyr::mutate(Pathways = factor(Pathways, levels = ptw.lvl)) %>%
  arrange(ptw.main.category, Pathways)

# all results (main pollutant + species + metals) ---------------------
tmp %<>% 
  dplyr::mutate(Windows = factor(Windows, levels = c("1y", "3y", "5y", "10y"))) %>% # set levels of windows
  dplyr::filter(adjust.p < 0.05) %>%  # only included pathways with FDR q<0.2
  dplyr::filter((n_significant_metabolites/n_metabolites) > 0.1) %>% # only included pathways with at least 10% significant metabolites
  dplyr::mutate(Expo_grp = str_replace_all(Expo_grp, "PM2.5_", "")) %>% 
  dplyr::mutate(Expo_grp = case_when(Expo_grp == "mixture_metal" ~ "Mixture(metal components)",
                                     Expo_grp == "mixture_pm25" ~ "Mixture(PM2.5 species)",
                                     Expo_grp == "PM25SO4" ~ "SO4",
                                     Expo_grp == "PM25NH4" ~ "NH4",
                                     Expo_grp == "PM25OC" ~ "OC",
                                     Expo_grp == "PM25EC" ~ "EC",
                                     TRUE ~ Expo_grp)) %>%
  dplyr::mutate(Expo_grp = factor(Expo_grp, levels = c("PM25", "O3", "OC", "EC", "SO4", "NH4", "Mixture(PM2.5 species)", 
                                                      "Bromine", "Copper", "Nickel", "Lead", "Potassium", "Silicon", "Vanadium", "Mixture(metal components)")))

# plot
p1 = tmp %>% 
  ggplot(aes(x = 0, y = Pathways)) + 
  geom_tile(aes(fill = ptw.main.category), alpha = 1) + 
  scale_fill_manual(values = ptw.color23[1:12]) + 
  scale_y_discrete(limits=rev, drop = TRUE) + 
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=16),
        legend.position = "right", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  guides(color= "none",  size="none")
p2 = tmp %>%
    ggplot()+
    geom_point(aes(x = Expo_grp, y = Pathways, color = Expo_grp, size = n_significant_metabolites), alpha = 1) + # +0.000001 to avoid p=0, which cause -log10(p) = Inf
    scale_size(range = c(3, 10)) + # set overall bubble size, which is different from size in above codeline
    facet_grid(. ~ Windows, scales = "free_x") + # set panel layout
    theme_bw() + 
    theme(axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x   = element_blank(),  # hide x labels
          axis.ticks.x  = element_blank(),
          strip.text = element_text(size = 16),
          legend.position = "right",
          ) + 
    scale_color_manual(values = expo.color22[1:15]) + 
    guides(color = guide_legend(override.aes = list(size = 10))) + 
    scale_y_discrete(limits=rev) + 
   # scale_color_manual(values = c("#d3d1cf", "#303233", "#4ab397", "#e9aece", "#4a9ecf", "#e0833b", "#ebe59c", "#903be0")) +
   #scale_x_continuous(limits = c(-1.5, 8.5)) + # add 1.5 to x-axis range to let bubbles be included in figure
    labs(title = "Pathway analysis results for PM2.5 metal components (FDR q<0.05) across time windows")
# merge
(p1|p2) + plot_layout(widths = c(0.2, 5), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right", 
        axis.text = element_text(size = 18))
ggsave(str_c(outpath, "/Visualizations/SEDAC/pathway analysis/raw_signal(p0.05)_as_input/Pathway_analysis_FDR0.05_ALL", ".png"), 
  height = 12, width = 20, dpi = 600) 


## 6-2:Pathway analysis viz- SEDAC (raw signal p<0.01 as input)-----------------------------------
file_names = list.files(path = str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.01"), pattern= "metapone_pathways_SEDAC_[a-zA-Z0-9._]+.xlsx$")
expo_names = file_names %>% str_extract("(?<=metapone_pathways_SEDAC_)[a-zA-Z0-9._]+(?=.xlsx)")
map2(file_names, expo_names, 
  ~readxl::read_xlsx(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.01/", .x)) %>% 
    dplyr::mutate(Exposure = .y) %>% 
    dplyr::mutate(Expo_grp = str_extract(Exposure, "[a-zA-Z0-9._]+(?=_\\d+y$)"),
                  Windows = str_extract(Exposure, "(?<=_)\\d+y$")) %>% 
    dplyr::select(Exposure, Expo_grp, Windows, everything())) %>% 
  rbindlist() -> plotdata

file_names = list.files(path = str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.01/TRAP_mixture"), pattern= "metapone_pathways_mixture_sedac.[a-zA-Z0-9._]+.xlsx$")
expo_names = file_names %>% str_extract("(?<=metapone_pathways_mixture_sedac.)[a-zA-Z0-9._]+(?=.xlsx)")
map2(file_names, expo_names, 
  ~readxl::read_xlsx(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.01/TRAP_mixture/", .x)) %>% 
    dplyr::mutate(Exposure = str_c("mixture_", .y)) %>% 
    dplyr::mutate(Expo_grp = str_extract(Exposure, "[a-zA-Z0-9._]+(?=_\\d+y$)"),
                  Windows = str_extract(.y, "(?<=_)\\w+$")) %>% 
    dplyr::select(Exposure, Expo_grp, Windows, everything())) %>% 
  rbindlist() -> plotdata2

plt.data = rbind(plotdata, plotdata2)
# plt.data %>% View()

# unique pathways, there were a total of 168 pathways identified among all group settings 
plt.data$Pathways %>% unique() %>% length()

###Plot----------------------------
#deal with p=0
plt.data$P <- ifelse(plt.data$p_value == 0, 0.0000001, plt.data$p_value)
plt.data %>% colnames()
# raw p<0.05
plt.data %>% 
  dplyr::mutate(P = if_else(p_value == 0, 0.0000001, p_value))  -> tmp

# export for manually checking pathway names
tmp %>% 
  fwrite(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.01", 
               "/0_pathway_analysis_results.csv"))

#####################################################################
# read modified tmp back to environment
tmp = fread(str_c(outpath, "/DataAnalysis_2_PathwayAnalysis_XtoM/Sensitive_analyses_using maximum datasets/raw0.01", 
               "/0_pathway_analysis_results_[manually modified].csv"))
#####################################################################
library(patchwork)

tmp %>% 
  dplyr::arrange(ptw.main.category) %>% 
  .$ptw.main.category %>% 
  unique() -> ptw.cate.lvl # pathway main category levels, n=12

tmp %>% 
  dplyr::select(Pathways, ptw.main.category) %>% 
  unique() %>% 
  arrange(ptw.main.category) %>% 
  .$Pathways %>% 
  unique()  -> ptw.lvl # pathway levels

ptw.color23 = c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
  "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
  "#393B79", "#637939", "#8C6D31"
)
expo.color22 = c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#97a884",
  "#E6AB02", "#A6761D", "#a19292", "#A6CEE3", "#535551",
  "#FDBF6F", "#CAB2D6"
)

# rank pathways by: 1) pathway main category; and 2) alphabet name order
colnames(tmp)
tmp %<>% 
  dplyr::mutate(Pathways = factor(Pathways, levels = ptw.lvl)) %>%
  arrange(ptw.main.category, Pathways)

# all results (main pollutant + species + metals) ---------------------
tmp %<>% 
  dplyr::mutate(Windows = factor(Windows, levels = c("1y", "3y", "5y", "10y"))) %>% # set levels of windows
  dplyr::filter(adjust.p < 0.05) %>%  # only included pathways with FDR q<0.2
  dplyr::filter((n_significant_metabolites/n_metabolites) > 0.1) %>% # only included pathways with at least 10% significant metabolites
  dplyr::mutate(Expo_grp = str_replace_all(Expo_grp, "PM2.5_", "")) %>% 
  dplyr::mutate(Expo_grp = case_when(Expo_grp == "mixture_metal" ~ "Mixture(metal components)",
                                     Expo_grp == "mixture_pm25" ~ "Mixture(PM2.5 species)",
                                     Expo_grp == "PM25SO4" ~ "SO4",
                                     Expo_grp == "PM25NH4" ~ "NH4",
                                     Expo_grp == "PM25OC" ~ "OC",
                                     Expo_grp == "PM25EC" ~ "EC",
                                     TRUE ~ Expo_grp)) %>%
  dplyr::mutate(Expo_grp = factor(Expo_grp, levels = c("PM25", "O3", "OC", "EC", "SO4", "NH4", "Mixture(PM2.5 species)", 
                                                      "Bromine", "Copper", "Nickel", "Lead", "Potassium", "Silicon", "Vanadium", "Mixture(metal components)")))

# plot
p1 = tmp %>% 
  ggplot(aes(x = 0, y = Pathways)) + 
  geom_tile(aes(fill = ptw.main.category), alpha = 1) + 
  scale_fill_manual(values = ptw.color23[1:12]) + 
  scale_y_discrete(limits=rev, drop = TRUE) + 
  theme_minimal() +
  theme(axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=16),
        legend.position = "right", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  guides(color= "none",  size="none")
p2 = tmp %>%
    ggplot()+
    geom_point(aes(x = Expo_grp, y = Pathways, color = Expo_grp, size = n_significant_metabolites), alpha = 1) + # +0.000001 to avoid p=0, which cause -log10(p) = Inf
    scale_size(range = c(3, 10)) + # set overall bubble size, which is different from size in above codeline
    facet_grid(. ~ Windows, scales = "free_x") + # set panel layout
    theme_bw() + 
    theme(axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x   = element_blank(),  # hide x labels
          axis.ticks.x  = element_blank(),
          strip.text = element_text(size = 16),
          legend.position = "right",
          ) + 
    scale_color_manual(values = expo.color22[1:15]) + 
    guides(color = guide_legend(override.aes = list(size = 10))) + 
    scale_y_discrete(limits=rev) + 
   # scale_color_manual(values = c("#d3d1cf", "#303233", "#4ab397", "#e9aece", "#4a9ecf", "#e0833b", "#ebe59c", "#903be0")) +
   #scale_x_continuous(limits = c(-1.5, 8.5)) + # add 1.5 to x-axis range to let bubbles be included in figure
    labs(title = "Pathway analysis results for PM2.5 metal components (FDR q<0.05) across time windows")
# merge
(p1|p2) + plot_layout(widths = c(0.2, 5), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right", 
        axis.text = element_text(size = 18))
ggsave(str_c(outpath, "/Visualizations/SEDAC/pathway analysis/raw_signal(p0.01)_as_input/Pathway_analysis_FDR0.05_ALL", ".png"), 
  height = 12, width = 20, dpi = 600) 


# 7- Sankey plot -------------------------------
library(ggsankey)

# load data -
metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))

metabolite[, c(1:3,6,7)] %>% 
  set_names(c("mz", "time", "mode", "metabolites", "exposures")) %>% 
  separate_rows(exposures, sep = ";\\s*") -> df

df %>% 
  dplyr::select(metabolites, exposures) %>% 
  dplyr::mutate(expo = str_extract(exposures, "\\w+(?= \\(\\d+y\\)$)"), 
                time = str_extract(exposures, "\\d+y")) -> df2
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Glucosamine",
              "Vitamin D2")
expo.lvl = c("Lead", "Silicon", "Vanadium", "O3", "SO4")
df2 %>% 
  dplyr::select(metabolites, expo) %>% 
  unique() %>% 
  dplyr::mutate(expo = factor(expo, levels = expo.lvl),
                metabolites = factor(metabolites, levels = meta.lvl)) %>% 
  make_long(expo, metabolites) -> df3

color6 = c("#71c07c", "#f0be54", "#6169dab6", "#ec8fc2", "#4ec517", "#999999")

df3 %>% 
  dplyr::mutate(node = factor(node, levels = c(expo.lvl, meta.lvl, NA)), # variable 'node' contains all levels across metabolites and pathways.
                next_node = factor(next_node, levels = c(meta.lvl, NA))) %>% # variable next_node only contains pathways of different levels.
  dplyr::mutate(fill_grp = case_when(
    x == "expo" ~ node,
    x == "metabolites" ~ "metabolites"
  )) %>% 
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(fill_grp), label = node)) +
  geom_sankey(alpha = 0.5) + 
  geom_sankey_label(size = 3.5, color = 1, fill = "white", 
                     hjust = 0.5, vjust = 0.5) +
  scale_fill_manual(values = color6) +
  theme_sankey() + 
  labs(title = "TRAP associated Level-1 metabolites (FDR q<0.2)")
ggsave(str_c(getwd(), "/0_Manuscript/figures/MWAS", "/Sankey_plot_with_lavel.png"), 
       height = 14, width = 8, dpi = 600)

df3 %>% 
  dplyr::mutate(node = factor(node, levels = c(expo.lvl, meta.lvl, NA)), # variable 'node' contains all levels across metabolites and pathways.
                next_node = factor(next_node, levels = c(meta.lvl, NA))) %>% # variable next_node only contains pathways of different levels.
  dplyr::mutate(fill_grp = case_when(
    x == "expo" ~ node,
    x == "metabolites" ~ "metabolites"
  )) %>% 
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(fill_grp), label = node)) +
  geom_sankey(alpha = 0.5) + 
  # geom_sankey_label(size = 3.5, color = 1, fill = "white", 
  #                    hjust = 0.5, vjust = 0.5) +
  scale_fill_manual(values = color6) +
  theme_sankey() + 
  labs(title = "TRAP associated Level-1 metabolites (FDR q<0.2)")
ggsave(str_c(getwd(), "/0_Manuscript/figures/MWAS", "/Sankey_plot_no_lavel.png"), 
       height = 14, width = 8, dpi = 600)


# 8. Sensitivity analysis ---------------------------------------------------------------
# import data 
purrr::map(list.files(glue("{outpath}/DataAnalysis_6_Sensitivity_analysis")), 
  ~ fread(glue("{outpath}/DataAnalysis_6_Sensitivity_analysis/{.x}"))) %>% 
  rbindlist() %>% # remove duplicate main analysis results
  unique() -> df

metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))
df
metabolite %<>% 
  dplyr::select(Mode, Components, "Brief name") %>% 
  unique() %>% 
  dplyr::rename(Metabolite = "Brief name")

df %>% 
  dplyr::select(Group, Exposure, Components, Mode, aCoef, Pvalue_a) %>% 
  left_join(metabolite, by = c("Components", "Mode")) %>% 
  dplyr::select(-Components, -Mode) %>% 
  dplyr::mutate(pollu = str_replace_all(Exposure, "PM2.5_|PM25", "")) %>% 
  dplyr::mutate(window = str_extract(pollu, "(?<=_)\\d+y$"), 
                pollu = str_remove(pollu, "_\\d+y$"), 
                pair = str_c(pollu, "_", window, "_", Metabolite)) %>% 
  dplyr::select(Group, aCoef, Pvalue_a, pollu, window, Metabolite, pair) -> plt.data
# unique levels
plt.data$window %>% unique()
plt.data$pollu %>% unique()
plt.data$Metabolite %>% unique()
# reset levels arranged by metabolites
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Mannosamine",
              "Vitamin D2")

plt.data %<>% 
  dplyr::mutate(Metabolite = factor(Metabolite, levels = meta.lvl), 
                window = factor(window, levels = plt.data$window %>% unique())) %>% 
  arrange(Metabolite, window, pollu)
plt.data$pair %>% unique() -> pair.lvl

plt.data %<>% 
  dplyr::mutate(pair = factor(pair, levels = pair.lvl))


# colors
meta.color = c(
  "#4A6FA5", "#C97A7A", "#6E7F5F", "#9A85C9", "#D9A066",
  "#8AB4D0", "#B85C57", "#9BAF7F", "#7765A8", "#E8BB7F",
  "#6B8FBF", "#8C3F3F", "#556B2F", "#C2AEDB", "#F0D9A8",
  "#A49E8D", "#69645C"
)
# p1: time windows
p1 = plt.data %>% 
  ggplot(aes(x = 0, y = pair)) + 
  geom_tile(aes(fill = window)) + 
  scale_fill_manual(values = c("#c0bdd1", "#8881BC", "#4850A0", "#103B95")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# p2: pollutants
p2 = plt.data %>% 
  ggplot(aes(x = 0, y = pair)) + 
  geom_tile(aes(fill = pollu)) + 
  scale_fill_manual(values = c("#e58080", "#6fcab4", "#7f5b9c", "#a5af4c", "#a1a19f")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# p3: metabolites
p3 = plt.data %>% 
  ggplot(aes(x = 0, y = pair)) + 
  geom_tile(aes(fill = Metabolite)) + 
  scale_fill_manual(values = meta.color) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# p4: main figure
p4 = plt.data %>%       
  ggplot(aes(Group, pair)) +
  geom_tile(aes(fill = aCoef)) + 
  scale_fill_gradient2(low = "#16a0ca", high = "#dd4444ed", mid = "white", midpoint = 0) + 
  geom_point(data = plt.data %>% dplyr::filter(Pvalue_a <0.05),  # remember change datasets here
             aes(Group, pair), 
             shape = 42, size = 8, color = "#4e4c4cf5") + 
  theme_bw() + 
  theme(panel.spacing = unit(0, "cm", data = NULL), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom",
        axis.title = element_text(size = 14), 
        plot.title = element_text(size = 14, face = "bold"))

# COMBINE
(p1|p2|p3|p4) + plot_layout(widths = c(0.5,0.5,0.5,5), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right", 
        axis.text = element_text(size = 18))
# export
ggsave(str_c(outpath, "/Visualizations/SEDAC/SM_Figxx", ".png"), 
       height = 15, width = 10, dpi = 600)


# 9- Stratification analysis and interaction analysis ----------------------------------------
library(patchwork)
library(RColorBrewer)
## 9.1- ApoE ----
fread(str_c(outpath, "/DataAnalysis_4_Stratification",
            "/SEDAC_individual_pollutant_stratified_by_ApoE_cat.csv")) %>% 
  dplyr::select(Strat, Exposure, Components, Mode, aCoef, LCI_a, UCI_a, Pvalue_a) %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"), 
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  arrange(Window, Exposure, Components, Mode) -> df1

# read metabolites information
metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))
metabolite %<>% 
  dplyr::select(Mode, Components, "Brief name") %>% 
  unique() %>% 
  dplyr::rename(Metabolite = "Brief name")
# combine
df1 %<>% 
  dplyr::left_join(metabolite, by = c("Components", "Mode"))

# pairs to be used
df1$Exposure %>% unique()
df1$Window %>% unique()
df1$Metabolite %>% unique()
# reset levels arranged by metabolites
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Mannosamine",
              "Vitamin D2")
# re-assign levels
df1 %>% 
  dplyr::mutate(Pair = str_c(Window, Exposure, Metabolite, sep = "_"), 
                Strat_grp = factor(Strat, levels = c("E4 absence", "E4 existance", "All")),
                Exposure = factor(Exposure, levels = c("O3", "SO4", "Silicon", "Lead", "Vanadium")), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y")), 
                Metabolite = factor(Metabolite, levels = meta.lvl)) %>% 
  arrange(Metabolite, Window, Exposure, Strat_grp) -> plt.dt.1
plt.dt.1$Pair %>% unique() -> pair.lvl
plt.dt.1 %<>% 
  dplyr::mutate(Pair = factor(Pair, levels = pair.lvl))

# import interaction analysis results
fread(str_c(outpath, "/DataAnalysis_5_Interaction",
            "/SEDAC_individual_pollutant_interaction_with_ApoE_cat.csv")) %>% 
  dplyr::filter(Group == "Interaction model-interaction effect") %>% 
  dplyr::filter(Pvalue_a <0.05) %>% 
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"),
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::select(Exposure, Window, Components, Mode) %>% 
  dplyr::mutate(flag="*", Strat_grp="All") -> inter
# combine interaction info
plt.dt.1 %<>% 
  left_join(inter, by = c("Exposure", "Window", "Components", "Mode", "Strat_grp"))

plt.dt.1 %<>% 
  dplyr::mutate(Strat_grp = factor(Strat_grp, levels = c("E4 absence", "E4 existance", "All")))

# sub-figures 
# colors
meta.color = c(
  "#4A6FA5", "#C97A7A", "#6E7F5F", "#9A85C9", "#D9A066",
  "#8AB4D0", "#B85C57", "#9BAF7F", "#7765A8", "#E8BB7F",
  "#6B8FBF", "#8C3F3F", "#556B2F", "#C2AEDB", "#F0D9A8",
  "#A49E8D", "#69645C"
)

p1 = plt.dt.1 %>% # time windows
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Window)) + 
  scale_fill_manual(values = c("#c0bdd1", "#8881BC", "#4850A0", "#103B95")) + 
  theme_minimal() +
  theme(#axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p2 = plt.dt.1 %>% # exposures
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Exposure)) +
  scale_fill_manual(values = c("#e58080", "#6fcab4", "#7f5b9c", "#a5af4c", "#a1a19f")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p3 = plt.dt.1 %>% 
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Metabolite)) +
  scale_fill_manual(values = meta.color) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p4 = plt.dt.1 %>%
  ggplot(aes(x = aCoef, y = Pair)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "#525050") +
  geom_point(aes(color = Strat_grp),
             size = 3,
             position = position_dodge(0.7)) + 
  geom_errorbar(aes(xmin = LCI_a, xmax = UCI_a, color = Strat_grp),
                width = 0.4,
                position = position_dodge(0.7)) + 
  scale_color_manual(values = c("#909e3d", "#4d5bac", "#46463e")) +
  # scale_x_continuous(limits = c(-6.5, 6), breaks = c(-5, -2.5, 0, 2.5, 5)) +
  geom_text(data = plt.dt.1,  # remember change datasets here
             aes(UCI_a+0.5, Pair, group = Strat_grp, label = flag), 
             size = 8, position = position_dodge(0.5), color = "red") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position = "right") + 
  facet_grid(.~Mode, scales = "free_x") + 
  guides(color = guide_legend(override.aes = list(size = 12)))

# merge
(p1|p2|p3|p4) + plot_layout(widths = c(0.5,0.5,0.5,7), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right",
        legend.key.size = unit(0.9, "cm"),
        axis.text = element_text(size = 18))
# save
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig4_ApoE_stratum", ".png"), 
       height = 20, width = 20, dpi = 600)


## 9.2- Age ----
fread(str_c(outpath, "/DataAnalysis_4_Stratification",
            "/SEDAC_individual_pollutant_stratified_by_Age_cat.csv")) %>% 
  dplyr::select(Strat, Exposure, Components, Mode, aCoef, LCI_a, UCI_a, Pvalue_a) %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"), 
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  arrange(Window, Exposure, Components, Mode) -> df2

# read metabolites information
metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))
metabolite %<>% 
  dplyr::select(Mode, Components, "Brief name") %>% 
  unique() %>% 
  dplyr::rename(Metabolite = "Brief name")
# combine
df2 %<>% 
  dplyr::left_join(metabolite, by = c("Components", "Mode"))

# pairs to be used
df2$Exposure %>% unique()
df2$Window %>% unique()
df2$Metabolite %>% unique()
# reset levels arranged by metabolites
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Mannosamine",
              "Vitamin D2")
# re-assign levels
df2 %>% 
  dplyr::mutate(Pair = str_c(Window, Exposure, Metabolite, sep = "_"), 
                Strat_grp = factor(Strat, levels = c("<70 yr", ">=70 yr", "All")),
                Exposure = factor(Exposure, levels = c("O3", "SO4", "Silicon", "Lead", "Vanadium")), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y")), 
                Metabolite = factor(Metabolite, levels = meta.lvl)) %>% 
  arrange(Metabolite, Window, Exposure, Strat_grp) -> plt.dt.2
plt.dt.2$Pair %>% unique() -> pair.lvl
plt.dt.2 %<>% 
  dplyr::mutate(Pair = factor(Pair, levels = pair.lvl))

# import interaction analysis results
fread(str_c(outpath, "/DataAnalysis_5_Interaction",
            "/SEDAC_individual_pollutant_interaction_with_Age_cat.csv")) %>% 
  dplyr::filter(Group == "Interaction model-interaction effect") %>% 
  dplyr::filter(Pvalue_a <0.05) %>% 
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"),
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::select(Exposure, Window, Components, Mode) %>% 
  dplyr::mutate(flag="*", Strat_grp="All") -> inter
# combine interaction info
plt.dt.2 %<>% 
  left_join(inter, by = c("Exposure", "Window", "Components", "Mode", "Strat_grp"))

plt.dt.2 %<>% 
  dplyr::mutate(Strat_grp = factor(Strat_grp, levels = c("<70 yr", ">=70 yr", "All")))

# sub-figures 
# colors
meta.color = c(
  "#4A6FA5", "#C97A7A", "#6E7F5F", "#9A85C9", "#D9A066",
  "#8AB4D0", "#B85C57", "#9BAF7F", "#7765A8", "#E8BB7F",
  "#6B8FBF", "#8C3F3F", "#556B2F", "#C2AEDB", "#F0D9A8",
  "#A49E8D", "#69645C"
)

p1 = plt.dt.2 %>% # time windows
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Window)) + 
  scale_fill_manual(values = c("#c0bdd1", "#8881BC", "#4850A0", "#103B95")) + 
  theme_minimal() +
  theme(#axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p2 = plt.dt.2 %>% # exposures
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Exposure)) +
  scale_fill_manual(values = c("#e58080", "#6fcab4", "#7f5b9c", "#a5af4c", "#a1a19f")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p3 = plt.dt.2 %>% 
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Metabolite)) +
  scale_fill_manual(values = meta.color) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p4 = plt.dt.2 %>%
  ggplot(aes(x = aCoef, y = Pair)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "#525050") +
  geom_point(aes(color = Strat_grp),
             size = 3,
             position = position_dodge(0.7)) + 
  geom_errorbar(aes(xmin = LCI_a, xmax = UCI_a, color = Strat_grp),
                width = 0.4,
                position = position_dodge(0.7)) + 
  scale_color_manual(values = c("#909e3d", "#4d5bac", "#46463e")) +
  # scale_x_continuous(limits = c(-6.5, 6), breaks = c(-5, -2.5, 0, 2.5, 5)) +
  geom_text(data = plt.dt.2,  # remember change datasets here --------------------
             aes(UCI_a+0.5, Pair, group = Strat_grp, label = flag), 
             size = 8, position = position_dodge(0.5), color = "red") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position = "right") + 
  facet_grid(.~Mode, scales = "free_x") + 
  guides(color = guide_legend(override.aes = list(size = 12)))

# merge
(p1|p2|p3|p4) + plot_layout(widths = c(0.5,0.5,0.5,7), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right",
        legend.key.size = unit(0.9, "cm"),
        axis.text = element_text(size = 18))
# save
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig4_Age_stratum", ".png"), 
       height = 20, width = 20, dpi = 600)


## 9.3- ABC ----
fread(str_c(outpath, "/DataAnalysis_4_Stratification",
            "/SEDAC_individual_pollutant_stratified_by_ABC_cat.csv")) %>% 
  dplyr::select(Strat, Exposure, Components, Mode, aCoef, LCI_a, UCI_a, Pvalue_a) %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"), 
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  arrange(Window, Exposure, Components, Mode) -> df2

# read metabolites information
metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))
metabolite %<>% 
  dplyr::select(Mode, Components, "Brief name") %>% 
  unique() %>% 
  dplyr::rename(Metabolite = "Brief name")
# combine
df2 %<>% 
  dplyr::left_join(metabolite, by = c("Components", "Mode"))

# pairs to be used
df2$Exposure %>% unique()
df2$Window %>% unique()
df2$Metabolite %>% unique()
# reset levels arranged by metabolites
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Mannosamine",
              "Vitamin D2")
# re-assign levels
df2 %>% 
  dplyr::mutate(Pair = str_c(Window, Exposure, Metabolite, sep = "_"), 
                Strat_grp = factor(Strat, levels = c("Low", "High", "All")),
                Exposure = factor(Exposure, levels = c("O3", "SO4", "Silicon", "Lead", "Vanadium")), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y")), 
                Metabolite = factor(Metabolite, levels = meta.lvl)) %>% 
  arrange(Metabolite, Window, Exposure, Strat_grp) -> plt.dt.2
plt.dt.2$Pair %>% unique() -> pair.lvl
plt.dt.2 %<>% 
  dplyr::mutate(Pair = factor(Pair, levels = pair.lvl))

# import interaction analysis results
fread(str_c(outpath, "/DataAnalysis_5_Interaction",
            "/SEDAC_individual_pollutant_interaction_with_ABC_cat.csv")) %>% 
  dplyr::filter(Group == "Interaction model-interaction effect") %>% 
  dplyr::filter(Pvalue_a <0.05) %>% 
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"),
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::select(Exposure, Window, Components, Mode) %>% 
  dplyr::mutate(flag="*", Strat_grp="All") -> inter
# combine interaction info
plt.dt.2 %<>% 
  left_join(inter, by = c("Exposure", "Window", "Components", "Mode", "Strat_grp"))

plt.dt.2 %<>% 
  dplyr::mutate(Strat_grp = factor(Strat_grp, levels = c("Low", "High", "All")))

# sub-figures 
# colors
meta.color = c(
  "#4A6FA5", "#C97A7A", "#6E7F5F", "#9A85C9", "#D9A066",
  "#8AB4D0", "#B85C57", "#9BAF7F", "#7765A8", "#E8BB7F",
  "#6B8FBF", "#8C3F3F", "#556B2F", "#C2AEDB", "#F0D9A8",
  "#A49E8D", "#69645C"
)

p1 = plt.dt.2 %>% # time windows
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Window)) + 
  scale_fill_manual(values = c("#c0bdd1", "#8881BC", "#4850A0", "#103B95")) + 
  theme_minimal() +
  theme(#axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p2 = plt.dt.2 %>% # exposures
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Exposure)) +
  scale_fill_manual(values = c("#e58080", "#6fcab4", "#7f5b9c", "#a5af4c", "#a1a19f")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p3 = plt.dt.2 %>% 
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Metabolite)) +
  scale_fill_manual(values = meta.color) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p4 = plt.dt.2 %>%
  ggplot(aes(x = aCoef, y = Pair)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "#525050") +
  geom_point(aes(color = Strat_grp),
             size = 3,
             position = position_dodge(0.7)) + 
  geom_errorbar(aes(xmin = LCI_a, xmax = UCI_a, color = Strat_grp),
                width = 0.4,
                position = position_dodge(0.7)) + 
  scale_color_manual(values = c("#909e3d", "#4d5bac", "#46463e")) +
  # scale_x_continuous(limits = c(-6.5, 6), breaks = c(-5, -2.5, 0, 2.5, 5)) +
  geom_text(data = plt.dt.2,  # remember change datasets here --------------------
             aes(UCI_a+0.5, Pair, group = Strat_grp, label = flag), 
             size = 8, position = position_dodge(0.5), color = "red") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position = "right") + 
  facet_grid(.~Mode, scales = "free_x") + 
  guides(color = guide_legend(override.aes = list(size = 12)))

# merge
(p1|p2|p3|p4) + plot_layout(widths = c(0.5,0.5,0.5,7), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right",
        legend.key.size = unit(0.9, "cm"),
        axis.text = element_text(size = 18))
# save
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig4_ABC_stratum", ".png"), 
       height = 20, width = 20, dpi = 600)


## 9.4- Braak ----
fread(str_c(outpath, "/DataAnalysis_4_Stratification",
            "/SEDAC_individual_pollutant_stratified_by_Braak_cat.csv")) %>% 
  dplyr::select(Strat, Exposure, Components, Mode, aCoef, LCI_a, UCI_a, Pvalue_a) %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"), 
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  arrange(Window, Exposure, Components, Mode) -> df2

# read metabolites information
metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))
metabolite %<>% 
  dplyr::select(Mode, Components, "Brief name") %>% 
  unique() %>% 
  dplyr::rename(Metabolite = "Brief name")
# combine
df2 %<>% 
  dplyr::left_join(metabolite, by = c("Components", "Mode"))

# pairs to be used
df2$Exposure %>% unique()
df2$Window %>% unique()
df2$Metabolite %>% unique()
# reset levels arranged by metabolites
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Mannosamine",
              "Vitamin D2")
# re-assign levels
df2 %>% 
  dplyr::mutate(Pair = str_c(Window, Exposure, Metabolite, sep = "_"), 
                Strat_grp = factor(Strat, levels = c("Lower stage", "Higher stage", "All")),
                Exposure = factor(Exposure, levels = c("O3", "SO4", "Silicon", "Lead", "Vanadium")), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y")), 
                Metabolite = factor(Metabolite, levels = meta.lvl)) %>% 
  arrange(Metabolite, Window, Exposure, Strat_grp) -> plt.dt.2
plt.dt.2$Pair %>% unique() -> pair.lvl
plt.dt.2 %<>% 
  dplyr::mutate(Pair = factor(Pair, levels = pair.lvl))

# import interaction analysis results
fread(str_c(outpath, "/DataAnalysis_5_Interaction",
            "/SEDAC_individual_pollutant_interaction_with_Braak_cat.csv")) %>% 
  dplyr::filter(Group == "Interaction model-interaction effect") %>% 
  dplyr::filter(Pvalue_a <0.05) %>% 
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"),
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::select(Exposure, Window, Components, Mode) %>% 
  dplyr::mutate(flag="*", Strat_grp="All") -> inter
# combine interaction info
plt.dt.2 %<>% 
  left_join(inter, by = c("Exposure", "Window", "Components", "Mode", "Strat_grp"))

plt.dt.2 %<>% 
  dplyr::mutate(Strat_grp = factor(Strat_grp, levels = c("Lower stage", "Higher stage", "All")))


p1 = plt.dt.2 %>% # time windows
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Window)) + 
  scale_fill_manual(values = c("#c0bdd1", "#8881BC", "#4850A0", "#103B95")) + 
  theme_minimal() +
  theme(#axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p2 = plt.dt.2 %>% # exposures
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Exposure)) +
  scale_fill_manual(values = c("#e58080", "#6fcab4", "#7f5b9c", "#a5af4c", "#a1a19f")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p3 = plt.dt.2 %>% 
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Metabolite)) +
  scale_fill_manual(values = meta.color) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p4 = plt.dt.2 %>%
  ggplot(aes(x = aCoef, y = Pair)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "#525050") +
  geom_point(aes(color = Strat_grp),
             size = 3,
             position = position_dodge(0.7)) + 
  geom_errorbar(aes(xmin = LCI_a, xmax = UCI_a, color = Strat_grp),
                width = 0.4,
                position = position_dodge(0.7)) + 
  scale_color_manual(values = c("#909e3d", "#4d5bac", "#46463e")) +
  # scale_x_continuous(limits = c(-6.5, 6), breaks = c(-5, -2.5, 0, 2.5, 5)) +
  geom_text(data = plt.dt.2,  # remember change datasets here --------------------
             aes(UCI_a+0.5, Pair, group = Strat_grp, label = flag), 
             size = 8, position = position_dodge(0.5), color = "red") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position = "right") + 
  facet_grid(.~Mode, scales = "free_x") + 
  guides(color = guide_legend(override.aes = list(size = 12)))

# merge
(p1|p2|p3|p4) + plot_layout(widths = c(0.5,0.5,0.5,7), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right",
        legend.key.size = unit(0.9, "cm"),
        axis.text = element_text(size = 18))
# save
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig4_Braak_stratum", ".png"), 
       height = 20, width = 20, dpi = 600)


## 9.5- CERAD ----
fread(str_c(outpath, "/DataAnalysis_4_Stratification",
            "/SEDAC_individual_pollutant_stratified_by_CERAD_cat.csv")) %>% 
  dplyr::select(Strat, Exposure, Components, Mode, aCoef, LCI_a, UCI_a, Pvalue_a) %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"), 
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  arrange(Window, Exposure, Components, Mode) -> df2

# read metabolites information
metabolite = readxl::read_excel(str_c(getwd(), "/0_Manuscript/tables", 
        "/PN_Identified metabolites (Level-1 evidence).xlsx"))
metabolite %<>% 
  dplyr::select(Mode, Components, "Brief name") %>% 
  unique() %>% 
  dplyr::rename(Metabolite = "Brief name")
# combine
df2 %<>% 
  dplyr::left_join(metabolite, by = c("Components", "Mode"))

# pairs to be used
df2$Exposure %>% unique()
df2$Window %>% unique()
df2$Metabolite %>% unique()
# reset levels arranged by metabolites
meta.lvl =  c("Phthalic anhydride",
              "ADP", "cAMP", "Adenosine", "L-isoleucine", "Methionine", "CoA", 
              "LysoPC(16:0)", "cis-11-eicosenoic acid", "Hexadecanol", "Malate", 
              "Gluconate-6P", "UDP-GlcA", "UDP-GalNAc", "ADP-glucose", "D-Mannosamine",
              "Vitamin D2")
# re-assign levels
df2 %>% 
  dplyr::mutate(Pair = str_c(Window, Exposure, Metabolite, sep = "_"), 
                Strat_grp = factor(Strat, levels = c("Sparse", "Frequent", "All")),
                Exposure = factor(Exposure, levels = c("O3", "SO4", "Silicon", "Lead", "Vanadium")), 
                Window = factor(Window, levels = c("1y", "3y", "5y", "10y")), 
                Metabolite = factor(Metabolite, levels = meta.lvl)) %>% 
  arrange(Metabolite, Window, Exposure, Strat_grp) -> plt.dt.2
plt.dt.2$Pair %>% unique() -> pair.lvl
plt.dt.2 %<>% 
  dplyr::mutate(Pair = factor(Pair, levels = pair.lvl))

# import interaction analysis results
fread(str_c(outpath, "/DataAnalysis_5_Interaction",
            "/SEDAC_individual_pollutant_interaction_with_CERAD_cat.csv")) %>% 
  dplyr::filter(Group == "Interaction model-interaction effect") %>% 
  dplyr::filter(Pvalue_a <0.05) %>% 
  dplyr::mutate(Window = str_extract(Exposure, "\\d+y$"),
                Exposure = str_remove(Exposure, "_\\d+y$")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM2.5_", "")) %>% 
  dplyr::mutate(Exposure = str_replace_all(Exposure, "PM25", "")) %>%
  dplyr::select(Exposure, Window, Components, Mode) %>% 
  dplyr::mutate(flag="*", Strat_grp="All") -> inter
# combine interaction info
plt.dt.2 %<>% 
  left_join(inter, by = c("Exposure", "Window", "Components", "Mode", "Strat_grp"))

plt.dt.2 %<>% 
  dplyr::mutate(Strat_grp = factor(Strat_grp, levels = c("Sparse", "Frequent", "All")))


p1 = plt.dt.2 %>% # time windows
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Window)) + 
  scale_fill_manual(values = c("#c0bdd1", "#8881BC", "#4850A0", "#103B95")) + 
  theme_minimal() +
  theme(#axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p2 = plt.dt.2 %>% # exposures
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Exposure)) +
  scale_fill_manual(values = c("#e58080", "#6fcab4", "#7f5b9c", "#a5af4c", "#a1a19f")) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p3 = plt.dt.2 %>% 
  ggplot(aes(x = 0, y = Pair)) + 
  geom_tile(aes(fill = Metabolite)) +
  scale_fill_manual(values = meta.color) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text.y = element_blank(), # remove y2 axis facet labels
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 14)))
p4 = plt.dt.2 %>%
  ggplot(aes(x = aCoef, y = Pair)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "#525050") +
  geom_point(aes(color = Strat_grp),
             size = 3,
             position = position_dodge(0.7)) + 
  geom_errorbar(aes(xmin = LCI_a, xmax = UCI_a, color = Strat_grp),
                width = 0.4,
                position = position_dodge(0.7)) + 
  scale_color_manual(values = c("#909e3d", "#4d5bac", "#46463e")) +
  # scale_x_continuous(limits = c(-6.5, 6), breaks = c(-5, -2.5, 0, 2.5, 5)) +
  geom_text(data = plt.dt.2,  # remember change datasets here --------------------
             aes(UCI_a+0.5, Pair, group = Strat_grp, label = flag), 
             size = 8, position = position_dodge(0.5), color = "red") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position = "right") + 
  facet_grid(.~Mode, scales = "free_x") + 
  guides(color = guide_legend(override.aes = list(size = 12)))

# merge
(p1|p2|p3|p4) + plot_layout(widths = c(0.5,0.5,0.5,7), guides = "collect") &
  theme(plot.margin = unit(c(0,0,0,0), "cm"), 
        legend.position = "right",
        legend.key.size = unit(0.9, "cm"),
        axis.text = element_text(size = 18))
# save
ggsave(str_c(outpath, "/Visualizations/SEDAC/Fig4_CERAD_stratum", ".png"), 
       height = 20, width = 20, dpi = 600)


