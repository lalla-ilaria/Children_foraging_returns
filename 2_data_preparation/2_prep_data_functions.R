#functions to create data lists to pass to stancode

#first define average ages for traps and shells cause we're gonna use it elsewhere
mean_age_shells <- mean(real_data$shell_ppl$age)
mean_age_traps <- mean(real_data$trap_ppl$age)

##########################################################################
#PREPARE DATA - AGE ONLY
##########################################################################
#prepares lists to be passed to stan models. Can prepare either shell, trap, or all_trap data for respectively
make_list_data_age <- function(data = real_data, foraging_type ){
    #SHELLS
      d_shellppl <- data$shell_ppl
      d_shells <- data$shells
      
      #keep only foraging data
      d_shellppl <- d_shellppl %>% filter(data == "shells")
      
      #add index variables
      #index and sort all individuals so we can loop across them
      d_shellppl$index_id <- as.integer(as.factor(d_shellppl$anonymeID))
      d_shellppl <- d_shellppl[order(d_shellppl$index_id),]
      #add index for individuals in the foraging data
      for ( i in 1:nrow(d_shells)){
        d_shells$index_id[i] <- d_shellppl$index_id[which ( d_shellppl$anonymeID == d_shells$anonymeID[i])]
      }
    #TRAPS
      d_trapppl <- data$trap_ppl
      d_traps <- data$traps
      d_traps <- data$traps %>% filter(exposure > 0)
      
      #keep only foraging data
      d_trapppl <- d_trapppl %>% filter(data == "traps")
      
      #add index variables
      #index and sort all individuals so we can loop across them
      d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
      d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
      for ( i in 1:nrow(d_traps)){
        d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
      }
      
      #index id of best actor
      best_guy <- d_trapppl$index_id[which(d_trapppl$anonymeID == 13212)]
      
      
      #ALL TRAPS
      d_all_trapppl <- data$trap_ppl
      d_all_traps <- data$all_traps
      
      #keep only foraging data
      d_all_trapppl <- d_all_trapppl %>% filter(data == "traps")
      d_all_traps <- d_all_traps[which(d_all_traps$lenght_hour >= 1), ]
      
      #add index variables
      #index and sort all individuals so we can loop across them
      d_all_trapppl$index_id <- as.integer(as.factor(d_all_trapppl$anonymeID))
      d_all_trapppl <- d_all_trapppl[order(d_all_trapppl$index_id),]
      for ( i in 1:nrow(d_all_traps)){
        d_all_traps$index_id[i] <- d_all_trapppl$index_id[which ( d_all_trapppl$anonymeID == d_all_traps$anonymeID[i])]
      }
      
      #index id of best actor
      best_all_trap_guy <- d_all_trapppl$index_id[which(d_all_trapppl$anonymeID == 13212)]
      
      
      dat_shells_age <- list(
        N = nrow(d_shellppl),
        M = nrow(d_shells),
        age = d_shellppl$age / mean(data$shell_ppl$age),
        returns = as.numeric(d_shells$returns)/1000,
        duration = d_shells$lenght_min/mean(d_shells$lenght_min),
        tide = d_shells$tide_avg_depth,
        ID_i= d_shells$index_id
      )
      
      
      dat_traps_age <- list(
        N = nrow(d_trapppl),                       #n individuals in total sample
        M = nrow(d_traps),                         #n trip/person
        ID_i= d_traps$index_id,                    #index of person of trip 
        success = d_traps$success,                 #whether trap captured something
        age = d_trapppl$age / mean(data$trap_ppl$age),
        duration = d_traps$exposure/mean(d_traps$exposure),
        best_guy = best_guy
      )
      
      dat_all_traps_age <- list(
        N = nrow(d_all_trapppl),                       #n individuals in total sample
        M = nrow(d_all_traps),                         #n trip/person
        ID_i= d_all_traps$index_id,                    #index of person of trip 
        success = d_all_traps$success,                 #whether trap captured something
        age = d_all_trapppl$age / mean(data$trap_ppl$age),
        duration = d_all_traps$lenght_hour/mean(d_all_traps$lenght_hour),
        best_all_trap_guy = best_all_trap_guy
      )
      
      #return either type of data  
  if(foraging_type == "shells"){
    return(dat_shells_age)
  }else{
    if(foraging_type == "traps"){
      return(dat_traps_age)
    } else{
      if(foraging_type == "all_traps"){
        return(dat_all_traps_age)
      }
    }
  }
}



##########################################################################
#PREPARE DATA - ALL PREDICTORS
##########################################################################
make_list_data_all <- function(data = real_data, foraging_type ){

      #SHELLS
      d_shellppl <- data$shell_ppl
      d_shells <- data$shells
      d_shell_k <- data$shell_k
      
      #add index variables
      #index and sort all individuals so we can loop across them
      d_shellppl$index_id <- as.integer(as.factor(d_shellppl$anonymeID))
      d_shellppl <- d_shellppl[order(d_shellppl$index_id),]
      #add index for individuals in the foraging data
      for ( i in 1:nrow(d_shells)){
        d_shells$index_id[i] <- d_shellppl$index_id[which ( d_shellppl$anonymeID == d_shells$anonymeID[i])]
      }
      #sort knowledge data
      for (i in 1:3) d_shell_k[[i]] <-  d_shell_k[[i]][ order(row.names(d_shell_k[[i]])), ]
      
      #TRAPS
      d_trapppl <- data$trap_ppl
      d_traps <- data$traps
      d_traps <- data$traps %>% filter(exposure > 0)
      d_trap_k <- data$trap_k
      
      #add index variables
      #index and sort all individuals so we can loop across them
      d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
      d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
      for ( i in 1:nrow(d_traps)){
        d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
      }
      
      #sort knowledge data
      for (i in 1:3) d_trap_k[[i]] <-  d_trap_k[[i]][ order(row.names(d_trap_k[[i]])), ]
      
      #ALL TRAPS
      d_all_trapppl <- data$trap_ppl
      d_all_traps <- data$all_traps
      d_all_trap_k <- data$trap_k
      
      #add index variables
      #index and sort all individuals so we can loop across them
      d_all_trapppl$index_id <- as.integer(as.factor(d_all_trapppl$anonymeID))
      d_all_trapppl <- d_all_trapppl[order(d_all_trapppl$index_id),]
      for ( i in 1:nrow(d_all_traps)){
        d_all_traps$index_id[i] <- d_all_trapppl$index_id[which ( d_all_trapppl$anonymeID == d_all_traps$anonymeID[i])]
      }
      
      #remove traps shorter than one hour
      d_all_traps <- d_all_traps[which(d_all_traps$lenght_hour >= 1), ]
      
      #sort knowledge data
      d_all_trap_k <- for (i in 1:3) d_trap_k[[i]][ order(row.names(d_trap_k[[i]])), ]
      
      #SHELLS
      dat_shells_all <- list(
        #foraging data
        N = nrow(d_shellppl),                       #n individuals in total sample
        M = nrow(d_shells),                         #n trip/person
        ID_i= d_shells$index_id,                    #index of person of trip 
        returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
        age = d_shellppl$age / mean(real_data$shell_ppl$age),
        sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
        duration = d_shells$lenght_min/mean(d_shells$lenght_min),
        tide = d_shells$tide_avg_depth,
        #height data
        has_height = ifelse(is.na(d_shellppl$height), 0, 1),# #vector of 0/1 for whether height has to be imputed
        height = d_shellppl$height/mean(d_shellppl$height, na.rm = TRUE),
        min_height = 50/mean(d_shellppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
        #grip data
        has_grip = ifelse(is.na(d_shellppl$grip), 0, 1),# #vector of 0/1 for whether grip has to be imputed
        grip = d_shellppl$grip/mean(d_shellppl$grip, na.rm = TRUE),
        #knowledge data
        has_knowledge = ifelse(is.na(d_shellppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
        knowledge_nit = d_shellppl$knowledge/mean(d_shellppl$knowledge, na.rm = TRUE),
        F_n = ncol(d_shell_k[[1]]),                        #n items in freelist
        Q_n = ncol(d_shell_k[[2]]),                        #n items in questionnaire
        R_n = ncol(d_shell_k[[3]]),                        #n items in image recognition
        answers_f = d_shell_k[[1]],                         #all answers from freelist
        answers_q = d_shell_k[[2]],                         #all answers from freelist
        answers_r = d_shell_k[[3]]                         #all answers from freelist
      )
      
      #TRAPS
      dat_traps_all <- list(
        #foraging data
        N = nrow(d_trapppl),                       #n individuals in total sample
        M = nrow(d_traps),                         #n trip/person
        ID_i= d_traps$index_id,                    #index of person of trip 
        has_foraging = ifelse(d_trapppl$data == "traps", 1, 0),
        success = d_traps$success,                 #whether trap captured something
        age = d_trapppl$age / mean(real_data$trap_ppl$age),
        sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
        duration = d_traps$exposure/mean(d_traps$exposure),
        #height data
        has_height = ifelse(is.na(d_trapppl$height), 0, 1),# #vector of 0/1 for whether height has to be imputed
        height = d_trapppl$height/mean(d_trapppl$height, na.rm = TRUE),
        min_height = 50/mean(d_trapppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
        #grip data
        has_grip = ifelse(is.na(d_trapppl$grip), 0, 1),# #vector of 0/1 for whether grip has to be imputed
        grip = d_trapppl$grip/mean(d_trapppl$grip, na.rm = TRUE),
        #knowledge data
        has_knowledge = ifelse(is.na(d_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
        knowledge_nit = d_trapppl$knowledge/mean(d_trapppl$knowledge, na.rm = TRUE),
        F_n = ncol(d_trap_k[[1]]),                        #n items in freelist
        Q_n = ncol(d_trap_k[[2]]),                        #n items in questionnaire
        R_n = ncol(d_trap_k[[3]]),                        #n items in image recognition
        answers_f = d_trap_k[[1]],                         #all answers from freelist
        answers_q = d_trap_k[[2]],                         #all answers from freelist
        answers_r = d_trap_k[[3]]                         #all answers from freelist
      )
      
      #ALL TRAPS
      dat_all_traps_all <- list(
        #foraging data
        N = nrow(d_all_trapppl),                       #n individuals in total sample
        M = nrow(d_all_traps),                         #n trip/person
        ID_i= d_all_traps$index_id,                    #index of person of trip 
        has_foraging = ifelse(d_all_trapppl$data == "traps", 1, 0),
        success = d_all_traps$success,                 #whether trap captured something
        age = d_all_trapppl$age / mean(data$trap_ppl$age),
        sex = ifelse(d_all_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
        duration = d_all_traps$lenght_hour/mean(d_all_traps$lenght_hour),
        #height data
        has_height = ifelse(is.na(d_all_trapppl$height), 0, 1),# #vector of 0/1 for whether height has to be imputed
        height = d_all_trapppl$height/mean(d_all_trapppl$height, na.rm = TRUE),
        min_height = 50/mean(d_all_trapppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
        #grip data
        has_grip = ifelse(is.na(d_all_trapppl$grip), 0, 1),# #vector of 0/1 for whether grip has to be imputed
        grip = d_all_trapppl$grip/mean(d_all_trapppl$grip, na.rm = TRUE),
        #knowledge data
        has_knowledge = ifelse(is.na(d_all_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
        knowledge_nit = d_all_trapppl$knowledge/mean(d_all_trapppl$knowledge, na.rm = TRUE),
        F_n = ncol(d_trap_k[[1]]),                        #n items in freelist
        Q_n = ncol(d_trap_k[[2]]),                        #n items in questionnaire
        R_n = ncol(d_trap_k[[3]]),                        #n items in image recognition
        answers_f = d_trap_k[[1]],                         #all answers from freelist
        answers_q = d_trap_k[[2]],                         #all answers from freelist
        answers_r = d_trap_k[[3]]                         #all answers from freelist
      )
      
      #return either type of data  
  if(foraging_type == "shells"){
    return(dat_shells_all)
  }else{
    if(foraging_type == "traps"){
      return(dat_traps_all)
    } else{
      if(foraging_type == "all_traps"){
        return(dat_all_traps_all)
      }
    }
  }
}

##########################################################################
#PREPARE DATA - COMPLETE CASES ONLY
##########################################################################
#prepares lists to be passed to stan models. Can prepare either shell, trap, or all_trap data for respectively
make_list_data_complete_cases <- function(data = real_data, foraging_type ){
  #SHELLS
  d_shellppl <- data$shell_ppl
  d_shells <- data$shells
  d_shell_k <- data$shell_k
  
  #keep only foraging data
  d_shellppl <- d_shellppl %>% filter(data == "shells")
  d_shellppl <- d_shellppl[complete.cases(d_shellppl),]
  
  #add index variables
  #index and sort all individuals so we can loop across them
  d_shellppl$index_id <- as.integer(as.factor(d_shellppl$anonymeID))
  d_shellppl <- d_shellppl[order(d_shellppl$index_id),]
  #add index for individuals in the foraging data
  d_shells <- d_shells[ which(d_shells$anonymeID %in% d_shellppl$anonymeID), ]
  for ( i in 1:nrow(d_shells)){
    d_shells$index_id[i] <- d_shellppl$index_id[which ( d_shellppl$anonymeID == d_shells$anonymeID[i])]
  }
  
  #filter and sort knowledge data
  for (i in 1:3) d_shell_k[[i]] <-  d_shell_k[[i]][ which(rownames(d_shell_k[[i]]) %in% d_shellppl$anonymeID), ]
  for (i in 1:3) d_shell_k[[i]] <-  d_shell_k[[i]][ order(row.names(d_shell_k[[i]])), ]
  
  #TRAPS
  d_trapppl <- data$trap_ppl
  d_traps <- data$traps
  d_traps <- data$traps %>% filter(exposure > 0)
  d_trap_k <- data$trap_k
  
  #keep only foraging data
  d_trapppl <- d_trapppl %>% filter(data == "traps")
  d_trapppl <- d_trapppl[complete.cases(d_trapppl),]
  
  #add index variables
  #index and sort all individuals so we can loop across them
  d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
  d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
  d_traps <- d_traps[ which(d_traps$anonymeID %in% d_trapppl$anonymeID), ]
  for ( i in 1:nrow(d_traps)){
    d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
  }
  
  for (i in 1:3) d_trap_k[[i]] <-  d_trap_k[[i]][ which(rownames(d_trap_k[[i]]) %in% d_trapppl$anonymeID), ]
  for (i in 1:3) d_trap_k[[i]] <-  d_trap_k[[i]][ order(row.names(d_trap_k[[i]])), ]
  
  
  #index id of best actor
  best_guy <- d_trapppl$index_id[which(d_trapppl$anonymeID == 13212)]
  
  
  #ALL TRAPS
  d_all_trapppl <- data$trap_ppl
  d_all_traps <- data$all_traps
  d_all_trap_k <- data$trap_k
  
  #keep only foraging data
  d_all_trapppl <- d_all_trapppl %>% filter(data == "traps")
  d_all_trapppl <- d_all_trapppl[complete.cases(d_all_trapppl),]
  d_all_traps <- d_all_traps[which(d_all_traps$lenght_hour >= 1), ]
  
  #add index variables
  #index and sort all individuals so we can loop across them
  d_all_trapppl$index_id <- as.integer(as.factor(d_all_trapppl$anonymeID))
  d_all_trapppl <- d_all_trapppl[order(d_all_trapppl$index_id),]
  d_all_traps <- d_all_traps[ which(d_all_traps$anonymeID %in% d_all_trapppl$anonymeID), ]
  for ( i in 1:nrow(d_all_traps)){
    d_all_traps$index_id[i] <- d_all_trapppl$index_id[which ( d_all_trapppl$anonymeID == d_all_traps$anonymeID[i])]
  }
  
  for (i in 1:3) d_trap_k[[i]] <-  d_trap_k[[i]][ which(rownames(d_trap_k[[i]]) %in% d_trapppl$anonymeID), ]
  for (i in 1:3) d_trap_k[[i]] <-  d_trap_k[[i]][ order(row.names(d_trap_k[[i]])), ]
  
  #index id of best actor
  best_all_trap_guy <- d_all_trapppl$index_id[which(d_all_trapppl$anonymeID == 13212)]
  
  
  dat_shells_cc <- list(
    N = nrow(d_shellppl),
    M = nrow(d_shells),
    age = d_shellppl$age / mean(data$shell_ppl$age),
    returns = as.numeric(d_shells$returns)/1000,
    duration = d_shells$lenght_min/mean(d_shells$lenght_min),
    tide = d_shells$tide_avg_depth,
    ID_i= d_shells$index_id,
    height = d_shellppl$height/mean(d_shellppl$height, na.rm = TRUE),
    min_height = 50/mean(d_shellppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
    #grip data
    grip = d_shellppl$grip/mean(d_shellppl$grip, na.rm = TRUE),
    #knowledge data
    knowledge_nit = d_shellppl$knowledge/mean(d_shellppl$knowledge, na.rm = TRUE),
    F_n = ncol(d_shell_k[[1]]),                        #n items in freelist
    Q_n = ncol(d_shell_k[[2]]),                        #n items in questionnaire
    R_n = ncol(d_shell_k[[3]]),                        #n items in image recognition
    answers_f = d_shell_k[[1]],                         #all answers from freelist
    answers_q = d_shell_k[[2]],                         #all answers from freelist
    answers_r = d_shell_k[[3]]                         #all answers from freelist
  )
  
  
  dat_traps_cc <- list(
    N = nrow(d_trapppl),                       #n individuals in total sample
    M = nrow(d_traps),                         #n trip/person
    ID_i= d_traps$index_id,                    #index of person of trip 
    success = d_traps$success,                 #whether trap captured something
    age = d_trapppl$age / mean(data$trap_ppl$age),
    duration = d_traps$exposure/mean(d_traps$exposure),
    best_guy = best_guy,
    height = d_trapppl$height/mean(d_trapppl$height, na.rm = TRUE),
    min_height = 50/mean(d_trapppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
    #grip data
    grip = d_trapppl$grip/mean(d_trapppl$grip, na.rm = TRUE),
    #knowledge data
    knowledge_nit = d_trapppl$knowledge/mean(d_trapppl$knowledge, na.rm = TRUE),
    F_n = ncol(d_trap_k[[1]]),                        #n items in freelist
    Q_n = ncol(d_trap_k[[2]]),                        #n items in questionnaire
    R_n = ncol(d_trap_k[[3]]),                        #n items in image recognition
    answers_f = d_trap_k[[1]],                         #all answers from freelist
    answers_q = d_trap_k[[2]],                         #all answers from freelist
    answers_r = d_trap_k[[3]]                         #all answers from freelist
  )
  
  dat_all_traps_cc <- list(
    N = nrow(d_all_trapppl),                       #n individuals in total sample
    M = nrow(d_all_traps),                         #n trip/person
    ID_i= d_all_traps$index_id,                    #index of person of trip 
    success = d_all_traps$success,                 #whether trap captured something
    age = d_all_trapppl$age / mean(data$trap_ppl$age),
    duration = d_all_traps$lenght_hour/mean(d_all_traps$lenght_hour),
    best_all_trap_guy = best_all_trap_guy,
    height = d_all_trapppl$height/mean(d_all_trapppl$height, na.rm = TRUE),
    min_height = 50/mean(d_all_trapppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
    #grip data
    grip = d_all_trapppl$grip/mean(d_all_trapppl$grip, na.rm = TRUE),
    #knowledge data
    knowledge_nit = d_all_trapppl$knowledge/mean(d_all_trapppl$knowledge, na.rm = TRUE),
    F_n = ncol(d_trap_k[[1]]),                        #n items in freelist
    Q_n = ncol(d_trap_k[[2]]),                        #n items in questionnaire
    R_n = ncol(d_trap_k[[3]]),                        #n items in image recognition
    answers_f = d_trap_k[[1]],                         #all answers from freelist
    answers_q = d_trap_k[[2]],                         #all answers from freelist
    answers_r = d_trap_k[[3]]                         #all answers from freelist
  )
  
  #return either type of data  
  if(foraging_type == "shells"){
    return(dat_shells_cc)
  }else{
    if(foraging_type == "traps"){
      return(dat_traps_cc)
    } else{
      if(foraging_type == "all_traps"){
        return(dat_all_traps_cc)
      }
    }
  }
}

