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
  dplyr::select(Case_Number, OBJECTID, GUID, date_of_death, POINT_X, POINT_Y) %>% 
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
fread(file = str_c(inpath, "/Air_pollution_data/CACES_LUR model.csv")) -> ap.data
ap.data$year %>% range(na.rm = T) # 1990 2020

# make lon/lat as sf objects
ap.data_sf <- st_as_sf(ap.data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
demo_sf <- st_as_sf(demo.ids, coords = c("POINT_X", "POINT_Y"), crs = 4326, remove = FALSE)


# 3- Assign air pollution data into demographic tables ----
# 3.1 function to extract nearest point
Ex.near.point = function(pollu){
  pollu_sub <- ap.data_sf %>% filter(pollutant == pollu)

  # unique locations, keyed on lon/lat only
  locs <- pollu_sub %>% distinct(lon, lat, .keep_all = TRUE)

  # get nearest point for each participant
  idx <- st_nearest_feature(demo_sf, locs)
  # distance 
  dist_m <- st_distance(demo_sf, locs[idx, ], by_element = TRUE)
  
  matched <- tibble(
    globalid = demo.ids$GUID,
    Case_Number = demo.ids$Case_Number,
    CACES_lon  = locs$lon[idx],
    CACES_lat  = locs$lat[idx], 
    dist_km  = as.numeric(dist_m) / 1000
  )

  # pull ALL calendar years for the matched location, joining on lon + lat
  pollu_sub %<>%
    st_drop_geometry() %>%
    select(lon, lat, year, pred_wght)
  
  matched %>%
    left_join(pollu_sub, by = c("CACES_lon" = "lon", "CACES_lat" = "lat")) %>%
    mutate(pollutant = pollu) -> res

  return(res)
}
# extract
pm25_long <- Ex.near.point("pm25")
pm10_long <- Ex.near.point("pm10")

exposure_long <- bind_rows(pm25_long, pm10_long)


# 3.2 function to extract average air pollution level for each individual and pollutant
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
       exposure_long %>% 
         dplyr::filter(globalid == id & year == death_year) %>%
         dplyr::filter(pollutant == pollu) %>% 
         dplyr::select(pred_wght) %>% 
         set_names(str_c(pollu, "_", window)) -> res
       res %<>% cbind(tibble(row = i, globalid = id)) %>% slice(1)
       return(res)
    }# end F_mean 
  } else if (window != "1y"){ # for other scenarios, we need to judge if these years were all available
      F_mean = function(i){
        # participants' globalid and death year
        id = demo.ids$GUID[[i]]
        death_year = demo.ids$death_yr[[i]] - 1 # one year prior to death year, as discussed in group meeting
        exposure_long %>% 
          dplyr::filter(globalid == id) %>% 
          dplyr::filter(pollutant == pollu) -> tmp.data
        if((death_year - len_window) >= 2000){
          tmp.data %>% 
            dplyr::select(globalid, year, pred_wght) %>% 
            dplyr::filter(year %in% seq(death_year - len_window, death_year)) %>% 
            summarise(across(all_of("pred_wght"), ~mean(.x, na.rm = FALSE)))  %>% 
            set_names(str_c(pollu, "_", window)) -> res
        }else if ((death_year - len_window) < 2000) {
           NA %>% as_tibble() %>%
            set_names(str_c(pollu, "_", window)) -> res
        }
        res %<>% cbind(tibble(row = i, globalid = id)) %>% slice(1)
        return(res)
      }# end F_mean
    }# end else if
  # End create function
  # Execute
  res = map_dfr(1:nrow(demo.ids), F_mean) %>% 
    dplyr::select(-row, -globalid) # for all participants (1:nrows)
  return(res)
}



# pollutant list
pollu.list = c("pm25", "pm10")

## 3.1 - #1_year before death ----
tictoc::tic()
purrr::map2_dfc(
      .x = rep("1y", length(pollu.list)), 
      .y = pollu.list, 
      Ex.avg.ap) -> avg_1y
nrow(avg_1y) == nrow(demo.ids)
tictoc::toc() # ~10 seconds
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


# 4- Merge demographic data with air pollution data------------------------------------
exposure_long %>% 
  dplyr::select(-year, -pred_wght, -pollutant) %>% unique() %>% 
  as.data.table() %>% 
  dplyr::select(-globalid, -Case_Number) -> tmp

# merge
demo.ids %>%
  cbind(tmp) %>% 
  cbind(avg_1y) %>% 
  cbind(avg_3y) %>% 
  cbind(avg_5y) %>% 
  cbind(avg_10y) -> demo.ap

# 5- Export data------------------------------------
demo.ap %>% 
  fwrite(str_c(outpath, "/Cleaned_data_Air_pollution(CACES).csv"))

demo.ap$pm10_1y %>% hist()

demo.ap %>% 
  dplyr::filter(dist_km<5)
