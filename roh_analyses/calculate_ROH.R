library(tidyverse)

dir = ""
setwd(dir)

file_list = list()
for (i in 1:40){
  new_dir = paste0(dir, "ROH_Result_Trial_", i)
  setwd(new_dir)
  file = paste0(new_dir, "/", list.files(pattern = "*csv"))
  file_list[[i]] = file
}

roh_manip = function(list, trial_number){
  ### becasuse some files size 0, bind_rows cannot work correctly
  ### size elimination first
  roh_file_bigger0  = list[[trial_number]][sapply(list[[trial_number]], file.size) > 72] 
  roh_pres = data.frame()
  for (file in roh_file_bigger0){
    temp = read.csv(file, header = T)
    roh_pres = roh_pres %>% 
      bind_rows(., temp) %>% 
      ### use length in Morgan,not centiMorgan !!!!!  
      filter(lengthM > 0.04) %>% 
      mutate(Ind = gsub("_merged", "", iid)) }
  ### produce ROH present individuals NROH & SROH summary
  roh_pres_sum = roh_pres %>% 
    group_by(Ind) %>% 
    summarise(nroh_M = n(), sroh_M = sum(lengthM))
  ### list all ind id found in given directory > to find out which ind have zero ROH
  all_id = list[[i]] %>%
    as_tibble %>% 
    mutate(Ind = sub(".*/(cch[^_]+)_.*", "\\1", value)) %>%  
    select(Ind)
  roh_0 = data.frame(Ind = setdiff(unique(all_id$Ind), unique(roh_pres$Ind)), 
                     nroh_M = 0, sroh_M = 0) 
  ### merge roh present and zero roh inds 
  final_roh = roh_pres_sum %>% 
    bind_rows(., roh_0) %>% 
    mutate(trial = trial_number, 
           label = paste(Ind, trial, sep = ".")) 
  return(final_roh)
}

### produce data frame and rename them
for (i in 1:40){
  temp = roh_manip(file_list, i)
  assign(
    x = paste("roh_final", i, sep = "_"), 
    value = temp, 
    envir = .GlobalEnv
  )
  temp 
}


### store data frames' in list without changing their formats. 
data_frames <- mget(paste0("roh_final_", 1:40))


### bind all data frames
roh_final_all <- do.call(rbind, data_frames) %>% 
  mutate_at(c("trial"), as.character)
