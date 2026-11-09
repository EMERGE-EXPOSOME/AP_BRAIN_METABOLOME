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
# path
work.dir = str_c(getwd(), "")
inpath <- str_c(work.dir, "/1_Original_data")
outpath <- str_c(work.dir, "/3_Outputs")

# 1- Stratification analysis and interaction analysis ----------------------------------------
## 1.1- ApoE ----
path    <- glue("{work.dir}/0_Manuscript/1_rebuttal_NN/3_Supplementary tables_Oct14.xlsx")
sd_path <- glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz_metabolite_SD.xlsx")
 
# ---- 1. Read the two side-by-side blocks ------------------------------------
raw <- read_excel(path, sheet = "Table S6_Stratum_APOE", .name_repair = "minimal")
 
strat <- raw[, 1:11] %>%
  set_names(c("Strat_grp", "Strat", "Exposure", "Metabolite", "Mode", "n_obs",
              "beta", "lci", "uci", "t", "p")) %>%
  mutate(p = as.numeric(p),                                # p-values are stored as text in Excel
         pair_id = (row_number() - 1) %/% 3 + 1)          # 3 rows per pair
 
inter <- raw[, 13:20] %>%
  set_names(c("Interaction_term", "Group", "n_obs", "beta", "lci", "uci", "t", "p")) %>%
  mutate(p = as.numeric(p),
         pair_id = (row_number() - 1) %/% 3 + 1)
 
# ---- 1b. Standardize: divide beta and CI by the metabolite SD -----------------
# SD table has the same 102 rows, in the same order, as the results table
sd_tbl <- read_excel(sd_path)

strat <- strat %>%
  mutate(Mode = sd_tbl$Mode,
         SD   = sd_tbl$SD,
         across(c(beta, lci, uci), ~ .x / SD))   # SD change in metabolite per SD of exposure
 
# ---- 2. Duplicated pairs: for those detected in both negative and positive modes,keep the negative-mode feature ---------------------
# as the results across modes were very similar in terms of p-values and effect estimates
# full results could be still accessed in supplementary tables
keep_ids <- strat %>%
  filter(Strat == "All") %>%
  arrange(Mode != "Negative") %>%                 # negative mode first
  distinct(Exposure, Metabolite, .keep_all = TRUE) %>%
  pull(pair_id)
 
p_int <- inter %>%
  filter(Group == "Interaction model-interaction effect", pair_id %in% keep_ids) %>%
  select(pair_id, p_int = p)
 
# ---- 3. Labels, pathways, ordering -------------------------------------------
met_map <- tribble(
  ~Metabolite,                                ~met_lab,                 ~pathway,
  "D-GLUCOSAMINE",                            "Glucosamine",            "Nucleotide sugar biosynthesis",
  "ADENOSINE-5'-DIPHOSPHOGLUCOSE",            "ADP-glucose",            "Nucleotide sugar biosynthesis",
  "URIDINE-DIPHOSPHATE-N-ACETYLGLUCOSAMINE",  "UDP-GlcNAc",             "Nucleotide sugar biosynthesis",
  "URIDINE 5'-DIPHOSPHOGLUCURONIC ACID",      "UDP-GlcA",               "Nucleotide sugar biosynthesis",
  "6-PHOSPHOGLUCONIC ACID",                   "Gluconate-6P",           "Nucleotide sugar biosynthesis",
  "(S)-MALATE",                               "Malate",                 "Nucleotide sugar biosynthesis",
  "HEXADECANOL",                              "Hexadecanol",            "Lipid metabolism",
  "FA20:1(GONDOIC ACID)",                     "Gondoic acid",           "Lipid metabolism",
  "LYSOPC(16:0)",                             "LysoPC(16:0)",           "Lipid metabolism",
  "COA",                                      "Co-enzyme A",            "Amino acid metabolism",
  "N-ACETYL-DL-METHIONINE",                   "Methionine",             "Amino acid metabolism",
  "L-ISOLEUCINE",                             "Isoleucine",             "Amino acid metabolism",
  "ADENOSINE",                                "Adenosine",              "Purine metabolism",
  "ADENOSINE 3',5'-CYCLIC MONOPHOSPHATE",     "cAMP",                   "Purine metabolism",
  "ADENOSINE 3',5'-DIPHOSPHATE",              "ADP",                    "Purine metabolism",
  "VITAMIN D2",                               "Vitamin D2",             "Other",
  "PHTHALIC ANHYDRIDE",                       "Phthalic anhydride",     "Other"
) %>%
  mutate(met_lab = factor(met_lab, levels = met_lab))
 
path_levels <- c("Nucleotide sugar biosynthesis", "Lipid metabolism",
                 "Amino acid metabolism", "Purine metabolism", "Other")
 
# plotmath labels so subscripts render on any device
exp_map <- c(Lead = "Pb", O3 = "O[3]", Silicon = "Si",
             SO4 = "SO[4]^{2*'-'}", Vanadium = "V")
 
stopifnot(all(unique(strat$Metabolite) %in% met_map$Metabolite))   # no unmapped names
 
dat <- strat %>%
  filter(pair_id %in% keep_ids) %>%
  left_join(met_map, by = "Metabolite") %>%
  left_join(p_int, by = "pair_id") %>%
  separate_wider_delim(Exposure, "_", names = c("exp", "win"), cols_remove = FALSE) %>%
  mutate(win_n   = parse_number(win),
         win_lab = paste(win_n, "yr"),
         exp_lab = exp_map[exp],
         pathway = factor(str_wrap(pathway, 14), levels = str_wrap(path_levels, 14)),
         Strat   = factor(Strat, levels = c("All", "E4 absence", "E4 existance"),
                          labels = c("Main analysis", "APOE e4 absence", "APOE e4 existence")),
         sig     = factor(if_else(p < 0.05, "p < 0.05", "Not significant"),
                          levels = c("p < 0.05", "Not significant")))
 
# one row per pair, ordered by metabolite (met_map order), then window (top = 1)
rows <- dat %>%
  distinct(pair_id, met_lab, exp, exp_lab, win_n, win_lab, p_int) %>%
  arrange(met_lab, win_n, exp) %>%
  mutate(y = rev(row_number())) %>%                    # first row at the top
  mutate(block = cumsum(met_lab != lag(met_lab, default = first(met_lab))),
         shade = if_else(block %% 2 == 1, "#F2F2F2", "transparent"))
 
dat <- dat %>% left_join(select(rows, pair_id, y), by = "pair_id")
 
# metabolite label once per block, centred on its rows
met_lab_pos <- rows %>% group_by(met_lab) %>% summarise(y = mean(y), .groups = "drop")
 
# ---- 4. Shared pieces --------------------------------------------------------
red  <- "#C0392B"
grey <- "grey55"
gap  <- unit(2, "mm")
 
bg_layer <- list(
  geom_blank(data = rows, aes(y = y - 0.5), inherit.aes = FALSE),   # pin row range in every panel
  geom_blank(data = rows, aes(y = y + 0.5), inherit.aes = FALSE),
  geom_rect(data = filter(rows, shade != "transparent"),
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            fill = "#F2F2F2", inherit.aes = FALSE)
)
 
y_scale <- scale_y_continuous(expand = c(0, 0))
 
base_theme <- theme(
  panel.background = element_blank(),
  panel.grid       = element_blank(),
  panel.spacing.y  = gap,
  axis.text        = element_text(size = 11, color = "black"),
  axis.title       = element_text(size = 12, color = "black"),
  strip.background = element_blank(),
  strip.text.y     = element_blank()
)
 
# ---- 5. Panels ---------------------------------------------------------------
# left: text columns (Exposure | Window | Metabolite), each header is a strip
# so it sits level with the forest-panel strips
text_col <- function(data, label, hdr, parse = FALSE, size = 3.6, color = "grey25") {
  ggplot(mutate(data, hdr = hdr), aes(0.02, y, label = .data[[label]])) +
    bg_layer +
    geom_text(hjust = 0, parse = parse, size = size, color = color) +
    facet_grid(. ~ hdr) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    y_scale +
    base_theme +
    theme(strip.text.x = element_text(size = 12, face = "bold", hjust = 0,
                                      margin = margin(4, 0, 4, 2)),
          axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0))     # no white gaps between text columns
}
 
p_exp <- text_col(rows, "exp_lab", "Exposure", parse = TRUE)
p_win <- text_col(rows, "win_lab", "Window")
p_met <- text_col(met_lab_pos, "met_lab", "Metabolite", size = 3.8, color = "black")
 
# middle: forest plots, one column per stratum, shared x-axis
p_for <- ggplot(dat, aes(beta, y)) +
  bg_layer +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lci, xmax = uci, color = sig), height = 0, linewidth = 0.6) +
  geom_point(aes(color = sig, fill = sig), shape = 21, size = 2.4, stroke = 0.8) +
  facet_grid(. ~ Strat) +
  scale_color_manual(values = c("p < 0.05" = red, "Not significant" = grey), name = NULL) +
  scale_fill_manual(values  = c("p < 0.05" = red, "Not significant" = "white"), name = NULL) +
  y_scale +
  labs(x = "Standardized effect estimate (95% CI)", y = NULL) +
  base_theme +
  theme(panel.border    = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.spacing.x = unit(3, "mm"),
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        strip.text.x    = element_text(size = 12, face = "bold"))
 
# right: interaction p-values
p_pint <- ggplot(mutate(rows, hdr = "P for interaction"), aes(0.5, y)) +
  bg_layer +
  geom_text(aes(label = if_else(p_int < 0.001, "<0.001", sprintf("%.3f", p_int)),
                color = p_int < 0.05, fontface = if_else(p_int < 0.05, "bold", "plain")),
            size = 3.6, show.legend = FALSE) +
  facet_grid(. ~ hdr) +   # header as strip, level with forest strips
  scale_color_manual(values = c(`TRUE` = red, `FALSE` = "grey25")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  y_scale +
  base_theme +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        axis.ticks   = element_blank(),
        axis.text    = element_blank(),
        axis.title   = element_blank())
 
# ---- 6. Combine --------------------------------------------------------------
p_all <- (p_exp | p_win | p_met | p_for | p_pint) +
  plot_layout(widths = c(0.75, 0.75, 1.8, 9, 1.5), guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 11))

ggsave(glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz/strat_APOE.png"), p_all, width = 15, height = 9, dpi = 600)




## 1.2- Other stratified factors ----
# ===================== RUN ONCE =====================
path    <- glue("{work.dir}/0_Manuscript/1_rebuttal_NN/3_Supplementary tables_Oct14.xlsx")
sd_path <- glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz_metabolite_SD.xlsx")
out_dir <- glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz")

# one SD per feature x exposure sample; independent of the stratifying factor
sd_key <- read_excel(sd_path) %>% distinct(Exposure, Metabolite, Mode, SD)
stopifnot(!anyDuplicated(sd_key[c("Exposure", "Metabolite", "Mode")]))

# met_map, path_levels, exp_map, red, grey, gap, y_scale, base_theme: keep them as they are

# make_bg() replaces the global bg_layer, so each figure gets its own row shading
make_bg <- function(rows) list(
  geom_blank(data = rows, aes(y = y - 0.5), inherit.aes = FALSE),
  geom_blank(data = rows, aes(y = y + 0.5), inherit.aes = FALSE),
  geom_rect(data = filter(rows, shade != "transparent"),
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            fill = "#F2F2F2", inherit.aes = FALSE)
)

text_col <- function(data, label, hdr, bg, parse = FALSE, size = 3.6, color = "grey25") {
  ggplot(mutate(data, hdr = hdr), aes(0.02, y, label = .data[[label]])) +
    bg +
    geom_text(hjust = 0, parse = parse, size = size, color = color) +
    facet_grid(. ~ hdr) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    y_scale +
    base_theme +
    theme(strip.text.x = element_text(size = 12, face = "bold", hjust = 0,
                                      margin = margin(4, 0, 4, 2)),
          axis.ticks = element_blank(), axis.text = element_blank(),
          axis.title = element_blank(), plot.margin = margin(0, 0, 0, 0))
}

# ===================== PER FACTOR =====================
# strata: values exactly as in the sheet's Strat column (excluding "All")
# labels: how they should appear as panel titles
plot_strat <- function(sheet, strata, labels, out_file) {
  k <- length(strata) + 1                               # rows per pair: All + each stratum

  raw <- read_excel(path, sheet = sheet, .name_repair = "minimal")

  strat <- raw[, 1:11] %>%
    set_names(c("Strat_grp", "Strat", "Exposure", "Metabolite", "Mode", "n_obs",
                "beta", "lci", "uci", "t", "p")) %>%
    mutate(p = as.numeric(p), pair_id = (row_number() - 1) %/% k + 1)

  inter <- raw[, 13:20] %>%
    set_names(c("Interaction_term", "Group", "n_obs", "beta", "lci", "uci", "t", "p")) %>%
    mutate(p = as.numeric(p), pair_id = (row_number() - 1) %/% k + 1)

  stopifnot(all(strat$Strat == rep(c("All", strata), nrow(strat) / k)))   # catches wrong strata/order

  # standardize by matching on keys, not row position
  strat <- strat %>%
    left_join(sd_key, by = c("Exposure", "Metabolite", "Mode")) %>%
    mutate(across(c(beta, lci, uci), ~ .x / SD))
  stopifnot(!anyNA(strat$SD))                           # every row found its SD

  keep_ids <- strat %>%
    filter(Strat == "All") %>%
    arrange(Mode != "Negative") %>%
    distinct(Exposure, Metabolite, .keep_all = TRUE) %>%
    pull(pair_id)

  p_int <- inter %>%
    filter(Group == "Interaction model-interaction effect", pair_id %in% keep_ids) %>%
    select(pair_id, p_int = p)
  if (anyDuplicated(p_int$pair_id))
    stop("More than one interaction term per pair (factor has >2 levels): ",
         "supply one overall interaction p-value per pair.")

  stopifnot(all(unique(strat$Metabolite) %in% met_map$Metabolite))

  dat <- strat %>%
    filter(pair_id %in% keep_ids) %>%
    left_join(met_map, by = "Metabolite") %>%
    left_join(p_int, by = "pair_id") %>%
    separate_wider_delim(Exposure, "_", names = c("exp", "win"), cols_remove = FALSE) %>%
    mutate(win_n   = parse_number(win),
           win_lab = paste(win_n, "yr"),
           exp_lab = exp_map[exp],
           Strat   = factor(Strat, levels = c("All", strata),
                            labels = c("Main analysis", labels)),
           sig     = factor(if_else(p < 0.05, "p < 0.05", "Not significant"),
                            levels = c("p < 0.05", "Not significant")))

  rows <- dat %>%
    distinct(pair_id, met_lab, exp, exp_lab, win_n, win_lab, p_int) %>%
    arrange(met_lab, win_n, exp) %>%
    mutate(y = rev(row_number()),
           block = cumsum(met_lab != lag(met_lab, default = first(met_lab))),
           shade = if_else(block %% 2 == 1, "#F2F2F2", "transparent"))

  dat <- dat %>% left_join(select(rows, pair_id, y), by = "pair_id")
  met_lab_pos <- rows %>% group_by(met_lab) %>% summarise(y = mean(y), .groups = "drop")
  bg <- make_bg(rows)

  p_exp <- text_col(rows, "exp_lab", "Exposure", bg, parse = TRUE)
  p_win <- text_col(rows, "win_lab", "Window", bg)
  p_met <- text_col(met_lab_pos, "met_lab", "Metabolite", bg, size = 3.8, color = "black")

  p_for <- ggplot(dat, aes(beta, y)) +
    bg +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lci, xmax = uci, color = sig), height = 0, linewidth = 0.6) +
    geom_point(aes(color = sig, fill = sig), shape = 21, size = 2.4, stroke = 0.8) +
    facet_grid(. ~ Strat) +
    scale_color_manual(values = c("p < 0.05" = red, "Not significant" = grey), name = NULL) +
    scale_fill_manual(values  = c("p < 0.05" = red, "Not significant" = "white"), name = NULL) +
    y_scale +
    labs(x = "Standardized effect estimate (95% CI)", y = NULL) +
    base_theme +
    theme(panel.border    = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.spacing.x = unit(3, "mm"),
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          strip.text.x    = element_text(size = 12, face = "bold"))

  p_pint <- ggplot(mutate(rows, hdr = "P for interaction"), aes(0.5, y)) +
    bg +
    geom_text(aes(label = if_else(p_int < 0.001, "<0.001", sprintf("%.3f", p_int)),
                  color = p_int < 0.05, fontface = if_else(p_int < 0.05, "bold", "plain")),
              size = 3.6, show.legend = FALSE) +
    facet_grid(. ~ hdr) +
    scale_color_manual(values = c(`TRUE` = red, `FALSE` = "grey25")) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    y_scale +
    base_theme +
    theme(strip.text.x = element_text(size = 12, face = "bold"),
          axis.ticks = element_blank(), axis.text = element_blank(),
          axis.title = element_blank())

  p_all <- (p_exp | p_win | p_met | p_for | p_pint) +
    plot_layout(widths = c(0.75, 0.75, 1.8, 3 * k, 1.5), guides = "collect") &   # 3 width units per panel
    theme(legend.position = "bottom", legend.text = element_text(size = 11))

  ggsave(out_file, p_all, width = 6 + 3 * k, height = 9, dpi = 600)
  invisible(p_all)
}

# ===================== CALLS =====================
plot_strat(sheet  = "Table S7_Stratum_Sex",             #sheet name
           strata = c("Male", "Female"),                #exact values in Strat column
           labels = c("Male", "Female"),
           out_file = glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz/strat_sex.png"))

plot_strat(sheet  = "Table S8_Stratum_Age",            
           strata = c("<70 yr", ">=70 yr"),             
           labels = c("<70 yr", ">=70 yr"),
           out_file = glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz/strat_age.png"))

plot_strat(sheet  = "Table S9_Stratum_ABC",            
           strata = c("Low", "High"),             
           labels = c("Low", "High"),
           out_file = glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz/strat_ABC.png"))

plot_strat(sheet  = "Table S10_Stratum_CERAD",            
           strata = c("Sparse", "Frequent"),             
           labels = c("Sparse", "Frequent"),
           out_file = glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz/strat_CERAD.png"))

plot_strat(sheet  = "Table S11_Stratum_Braak",            
           strata = c("Lower stage", "Higher stage"),             
           labels = c("Lower stage", "Higher stage"),
           out_file = glue("{work.dir}/0_Manuscript/1_rebuttal_NN/viz/strat_Braak.png"))


