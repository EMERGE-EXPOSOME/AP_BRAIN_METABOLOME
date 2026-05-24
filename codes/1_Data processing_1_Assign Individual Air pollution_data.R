rm(list=(ls()))
# 0-Prepare steps =====================================
# load package ---------
pacman::p_load("readxl", 
               "purrr", 
               "lubridate", 
               "dplyr", 
               "magrittr",
               "stringr", 
               "rlist", 
               "tidyverse", 
               "sf", 
               "data.table",  
               "ddpcr")
getwd() # "C:/Users/MREN30/OneDrive - Emory/1_Ongoing_project/2_TRAP_ADRD"
# set work directory
work.dir = str_c(getwd(), "") # define any work directory you like
inpath <- str_c(work.dir, 
                "/1_Original_data")
outpath <- str_c(work.dir,
                "/3_Outputs")

# 1- Load demographic data and extract Geo grid ID and death date------------------------
fread(file = str_c(inpath, "/Demographics/Brain_bank_Marla_clean_432obs.csv")) -> demo
# View(demo)
demo %>% 
  dplyr::select(Case_Number, OBJECTID, GUID, date_of_death) %>% 
  dplyr::mutate(death_yr = as.numeric(str_extract(date_of_death, "\\d{4}(?=-\\d{2}-\\d{2})")),
                death_mt = as.numeric(str_extract(date_of_death, "(?<=\\d{4}-)\\d{2}(?=-\\d{2})"))) -> demo.ids
demo.ids$death_mt %>% hist()
demo.ids$death_yr %>% hist()

demo.ids %>% 
  dplyr::filter(is.na(death_yr))

## 1.1 Re-assign "death year" ----
# 1.2- Extract individuals who died after 2000 and before 2020------------------------
demo.ids$GUID %>% unique() %>% length() # total participants 425
demo.ids  %<>% 
  dplyr::filter(!is.na(death_yr)) %>% 
  arrange(death_yr) # delete rows with NA death_yr
demo.ids$death_yr %>% range(na.rm = T) # 1999 2020, ranges of death year
demo.ids %<>% 
  dplyr::filter(death_yr >= 2001 & death_yr <= 2020) # keep participants who died between 2001 and 2020, since we only have TRAP data between 2000 and 2019
demo.ids %>% nrow() # N = 346

# 2- Load aggregated air pollution data--------------
fread(file = str_c(inpath, "/Air_pollution_data/Brain_exposure_matched_2025Feb18_from_Zhihao.csv")) -> ap.data
ap.data$year %>% range(na.rm = T) # 2000 2019

# 3- Assign air pollution data into demographic tables ----
i = 1
# function to extract average air pollution level for each individual and pollutant
Ex.avg.ap = function(window, pollu){
  # (length of windows) -1
  len_window = as.numeric(str_extract(window, "\\d{1,2}"))-1
  print(str_c("Windows of ", len_window+1, " years before death for ", pollu, " begin !"))
  # Create function
  # for 1 year before death, we don't need to again judge if this specific year's data are available 
  if(window == "1y"){
    F_mean = function(i){
       # participants' globalid and death year
       id = demo.ids$GUID[[i]]
       death_year = demo.ids$death_yr[[i]] - 1 # one year prior to death year
       ap.data %>% 
         dplyr::filter(globalid == id & year == death_year) %>% 
         dplyr::select(all_of(pollu)) %>% 
         set_names(str_c(pollu, "_", window)) -> res
       res %<>% cbind(tibble(row = i))
       return(res)
    }# end F_mean 
  } else if (window != "1y"){ # for other scenarios, we need to judge if these years were all available
      F_mean = function(i){
        # participants' globalid and death year
        id = demo.ids$GUID[[i]]
        death_year = demo.ids$death_yr[[i]] - 1 # one year prior to death year, as discussed in group meeting
        ap.data %>% 
          dplyr::filter(globalid == id) -> tmp.data
        if((death_year - len_window) >= 2000){
          tmp.data %>% 
            dplyr::select(globalid, year, all_of(pollu)) %>% 
            dplyr::filter(year %in% seq(death_year - len_window, death_year)) %>% 
            summarise(across(all_of(pollu), ~mean(.x, na.rm = FALSE)))  %>% 
            set_names(str_c(pollu, "_", window)) -> res
        }else if ((death_year - len_window) < 2000) {
           NA %>% as_tibble() %>%
            set_names(str_c(pollu, "_", window)) -> res
        }
        res %<>% cbind(tibble(row = i))
        return(res)
      }# end F_mean
    }# end else if
  # End create function
  # Execute
  res = map_dfr(1:nrow(demo.ids), F_mean) %>% 
    dplyr::select(-row) # for all participants (1:nrows)
  return(res)
}

# pollutant list
pollu.list = ap.data %>% dplyr::select(Rline_CO:SEDAC_PM2.5_Zinc) %>% colnames()  # a total of 21 TRAP related variables

## 3.1 - #1_year before death ----
tictoc::tic()
purrr::map2_dfc(
      .x = rep("1y", length(pollu.list)), 
      .y = pollu.list, 
      Ex.avg.ap) -> avg_1y
nrow(avg_1y) == nrow(demo.ids)
tictoc::toc() # ~60 seconds
# 3.2- #3_years before death ----
purrr::map2_dfc(
      .x = rep("3y", length(pollu.list)), 
      .y = pollu.list, 
      Ex.avg.ap) -> avg_3y
nrow(avg_3y) == nrow(demo.ids)
# 3.3- #5_years before death ----
purrr::map2_dfc(
      .x = rep("5y", length(pollu.list)), 
      .y = pollu.list, 
      Ex.avg.ap) -> avg_5y
nrow(avg_5y) == nrow(demo.ids)
# 3.4- #10_years before death ----
purrr::map2_dfc(
      .x = rep("10y", length(pollu.list)), 
      .y = pollu.list, 
      Ex.avg.ap) -> avg_10y
nrow(avg_10y) == nrow(demo.ids)
# 3.5- #15_years before death ----
purrr::map2_dfc(
      .x = rep("15y", length(pollu.list)), 
      .y = pollu.list, 
      Ex.avg.ap) -> avg_15y
nrow(avg_15y) == nrow(demo.ids)

# 4- Merge demographic data with air pollution data------------------------------------ 
demo.ids %>% 
  cbind(avg_1y) %>% 
  cbind(avg_3y) %>% 
  cbind(avg_5y) %>% 
  cbind(avg_10y) %>% 
  cbind(avg_15y) -> demo.ap

# 5- Export data------------------------------------
demo.ap %>% 
  fwrite(str_c(outpath, "/Cleaned_data_Air_pollution.csv"))
