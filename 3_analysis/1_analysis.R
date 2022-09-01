library(rethinking)
library(dplyr)
library(rlist)
real_data <- list.load("2_data_preparation/processed_data.RData")


##########################################################################
#PREPARE DATA 
##########################################################################
#SHELLS
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells

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
d_trapppl <- real_data$trap_ppl
d_traps <- real_data$traps

#keep only foraging data
d_trapppl <- d_trapppl %>% filter(data == "traps")
d_traps <- d_traps[which(d_traps$lenght_hour >= 1), ]

#add index variables
#index and sort all individuals so we can loop across them
d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
for ( i in 1:nrow(d_traps)){
  d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
}


##################################################
#FIT MODELS AGE ONLY
##################################################
#NB data passed for all individuals in sample, add in models a ifelse not to estimate phi for people who are not interesting
#SHELLS
#create data frame
dat_shells <- list(
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  age = d_shellppl$age / mean(d_shellppl$age),
  returns = as.numeric(d_shells$returns)/1000,
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  ID_i= d_shells$index_id
)
m_shell_age <- cstan( file= "models/1_shell_age.stan" , data=dat_shells , chains=3, cores = 3, iter = 2000 )


dat_traps <- list(
  N = nrow(d_trapppl),                       #n individuals in total sample
  M = nrow(d_traps),                         #n trip/person
  ID_i= d_traps$index_id,                    #index of person of trip 
  success = d_traps$success,                 #whether trap captured something
  age = (d_trapppl$age / mean(d_trapppl$age)),
  duration = d_traps$lenght_hour/mean(d_traps$lenght_hour)
)
m_trap_age <- cstan( file= "models/1_trap_age.stan" , data=dat_traps , chains=3, cores = 3, iter = 2000 )

prec_shells<-precis(m_shells_age,dept=3,prob=0.95)
prec_traps<-precis(m_traps_age,dept=3,prob=0.95)
write.csv(prec_shells,"4_outcomes/model_fit/shells_prec_age.csv")
write.csv(prec_traps,"4_outcomes/model_fit/traps_prec_age.csv")
post_s<-extract.samples(m_shells_age)
post_t<-extract.samples(m_traps_age)
save(post_s, file = "4_outcomes/model_fit/post_s_age.rda")
save(post_t, file = "4_outcomes/model_fit/post_t_age.rda")



##########################################################################
#PREPARE DATA 
##########################################################################
#SHELLS
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells
d_shell_k <- real_data$shell_k

#add index variables
#index and sort all individuals so we can loop across them
d_shellppl$index_id <- as.integer(as.factor(d_shellppl$anonymeID))
d_shellppl <- d_shellppl[order(d_shellppl$index_id),]
#add index for individuals in the foraging data
for ( i in 1:nrow(d_shells)){
  d_shells$index_id[i] <- d_shellppl$index_id[which ( d_shellppl$anonymeID == d_shells$anonymeID[i])]
}
#sort knowledge data
d_shell_k <- d_shell_k[ order(row.names(d_shell_k)), ]

#TRAPS
d_trapppl <- real_data$trap_ppl
d_traps <- real_data$traps
d_trap_k <- real_data$trap_k

#add index variables
#index and sort all individuals so we can loop across them
d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
for ( i in 1:nrow(d_traps)){
  d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
}

#remove traps shorter than one hour
d_traps <- d_traps[which(d_traps$lenght_hour >= 1), ]

#sort knowledge data
d_trap_k <- d_trap_k[ order(row.names(d_trap_k)), ]


##################################################
#FIT MODELS all data
##################################################
#SHELLS
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),                       #n individuals in total sample
  M = nrow(d_shells),                         #n trip/person
  ID_i= d_shells$index_id,                    #index of person of trip 
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = d_shellppl$age / mean(d_shellppl$age),
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
  #knowledge = d_shellppl$knowledge/mean(d_shellppl$knowledge, na.rm = TRUE)
  Q = ncol(d_shell_k),                        #n items in freelist
  answers = d_shell_k                         #all answers from freelist
)
dat_shells[["answers"]][is.na(dat_shells[["answers"]])] <- -999
#dat_shells[["knowledge"]][is.na(dat_shells[["knowledge"]])] <- -999
dat_shells[["height"]][is.na(dat_shells[["height"]])] <- -999
dat_shells[["grip"]][is.na(dat_shells[["grip"]])] <- -999


#TRAPS
dat_traps <- list(
  #foraging data
  N = nrow(d_trapppl),                       #n individuals in total sample
  M = nrow(d_traps),                         #n trip/person
  ID_i= d_traps$index_id,                    #index of person of trip 
  has_foraging = ifelse(d_trapppl$data == "traps", 1, 0),
  success = d_traps$success,                 #whether trap captured something
  age = d_trapppl$age / mean(d_trapppl$age),
  sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_traps$lenght_hour/mean(d_traps$lenght_hour),
  #height data
  has_height = ifelse(is.na(d_trapppl$height), 0, 1),# #vector of 0/1 for whether height has to be imputed
  height = d_trapppl$height/mean(d_trapppl$height, na.rm = TRUE),
  min_height = 50/mean(d_trapppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
  #grip data
  has_grip = ifelse(is.na(d_trapppl$grip), 0, 1),# #vector of 0/1 for whether grip has to be imputed
  grip = d_trapppl$grip/mean(d_trapppl$grip, na.rm = TRUE),
  #knowledge data
  has_knowledge = ifelse(is.na(d_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  #knowledge = d_trapppl$knowledge/mean(d_trapppl$knowledge, na.rm = TRUE)
  Q = ncol(d_trap_k),                        #n items in freelist
  answers = d_trap_k                       #all answers from freelist
)

dat_traps[["answers"]][is.na(dat_traps[["answers"]])] <- -999
#dat_traps[["knowledge"]][is.na(dat_traps[["knowledge"]])] <- -999
dat_traps[["height"]][is.na(dat_traps[["height"]])] <- -999
dat_traps[["grip"]][is.na(dat_traps[["grip"]])] <- -999

# m_shells_all <- cstan( file= "models/7_shells_estimate_all.stan" , data=dat_shells , 
#                        chains=3, cores = 3, iter = 2000 )
# 
# m_traps_all <- cstan( file= "models/7_trap_estimate_all.stan" , data=dat_traps , 
#                            chains=3, cores = 3, iter = 2000 )
# 
# #test out of log -need to have same qualitative results
# m_oll_shells <- cstan( file= "models/8_shells_outoflog_all.stan" , data=dat_shells , 
#                        chains=3, cores = 3, iter = 2000 )
# 
# m_oll_traps <- cstan( file= "models/8_trap_outoflog_all.stan" , data=dat_traps , 
#                       chains=3, cores = 3, iter = 2000 )

# m_shells_nit <- cstan( file= "models/11_shells_nitemsknowledge.stan" , data=dat_shells , 
#                        chains=3, cores = 3, iter = 400 )
# 
# m_traps_nit <- cstan( file= "models/11_traps_nitemsknowledge.stan" , data=dat_traps , 
#                       chains=3, cores = 3, iter = 400 )

m_shells_neg <- cstan( file= "models/9_shells_estimate_all_negative_effects.stan" , data=dat_shells , 
                       chains=3, cores = 3, iter = 400 )

m_traps_neg <- cstan( file= "models/9_traps_estimate_all_negative_effects.stan" , data=dat_traps , 
                      chains=3, cores = 3, iter = 400 )

prec_shells<-precis(m_shells_neg,dept=3,prob=0.95)
prec_traps<-precis(m_traps_neg,dept=3,prob=0.95)
write.csv(prec_shells,"4_outcomes/model_fit/shells_prec.csv")
write.csv(prec_traps,"4_outcomes/model_fit/traps_prec.csv")
post_s<-extract.samples(m_shells_neg)
post_t<-extract.samples(m_traps_neg)
save(post_s, file = "4_outcomes/model_fit/post_s.rda")
save(post_t, file = "4_outcomes/model_fit/post_t.rda")

# m_shells_v <- cstan( file= "models/9v_shells_estimate_all_negative_effects.stan" , data=dat_shells , 
#                        chains=3, cores = 3, iter = 400 )
# 
# m_traps_v <- cstan( file= "models/9v_traps_estimate_all_negative_effects.stan" , data=dat_traps , 
#                       chains=3, cores = 3, iter = 400 )
# 
#m_shells_gripsex <- cstan( file= "models/10_shell_gripbyage.stan" , data=dat_shells , 
#                       chains=1, cores = 1, iter = 400 )
##########################################################################
#TRAPS AS POISSON
##########################################################################
# #TRAPS
# d_trapppl <- real_data$trap_ppl
# d_traps <- real_data$individual_traps
# 
# #keep only foraging data
# d_trapppl <- d_trapppl %>% filter(data == "traps")
# 
# #add index variables
# #index and sort all individuals so we can loop across them
# d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
# d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
# for ( i in 1:nrow(d_traps)){
#   d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
# }
# 
# dat_traps <- list(
#   N = nrow(d_trapppl),                       #n individuals in total sample
#   M = nrow(d_traps),                         #n trip/person
#   ID_i= d_traps$index_id,                    #index of person of trip 
#   success = d_traps$success,                 #whether trap captured something
#   age = (d_trapppl$age / mean(d_trapppl$age)),
#   duration = d_traps$exposure/mean(d_traps$exposure)
# )
# m_trap_age <- cstan( file= "models/1_trap_age_poisson.stan" , data=dat_traps , chains=3, cores = 3, iter = 500 )
# 
# 
# ####################
# #traps as poisson
# ####################
# 
# #TRAPS
# d_trapppl <- real_data$trap_ppl
# d_traps <- real_data$individual_traps
# d_trap_k <- real_data$trap_k
# 
# #add index variables
# #index and sort all individuals so we can loop across them
# d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
# d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
# for ( i in 1:nrow(d_traps)){
#   d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
# }
# 
# #sort knowledge data
# d_trap_k <- d_trap_k[ order(row.names(d_trap_k)), ]
# 
# #TRAPS
# dat_traps <- list(
#   #foraging data
#   N = nrow(d_trapppl),                       #n individuals in total sample
#   M = nrow(d_traps),                         #n trip/person
#   ID_i= d_traps$index_id,                    #index of person of trip 
#   has_foraging = ifelse(d_trapppl$data == "traps", 1, 0),
#   success = d_traps$success,                 #whether trap captured something
#   age = d_trapppl$age / mean(d_trapppl$age),
#   sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
#   duration = d_traps$exposure/mean(d_traps$exposure),
#   #height data
#   has_height = ifelse(is.na(d_trapppl$height), 0, 1),# #vector of 0/1 for whether height has to be imputed
#   height = d_trapppl$height/mean(d_trapppl$height, na.rm = TRUE),
#   min_height = 50/mean(d_trapppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
#   #grip data
#   has_grip = ifelse(is.na(d_trapppl$grip), 0, 1),# #vector of 0/1 for whether grip has to be imputed
#   grip = d_trapppl$grip/mean(d_trapppl$grip, na.rm = TRUE),
#   #knowledge data
#   has_knowledge = ifelse(is.na(d_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
#   Q = ncol(d_trap_k),                        #n items in freelist
#   answers = d_trap_k                       #all answers from freelist
# )
# 
# dat_traps[["answers"]][is.na(dat_traps[["answers"]])] <- -999
# dat_traps[["height"]][is.na(dat_traps[["height"]])] <- -999
# dat_traps[["grip"]][is.na(dat_traps[["grip"]])] <- -999
# 
# m_traps_pois <- cstan( file= "models/10_traps_poisson.stan" , data=dat_traps , 
#                       chains=3, cores = 3, iter = 2000 )

##################
#TIDE MODEL
##################

for(i in 1:nrow(d_shells)){
  d_shells$age[i] <- d_shellppl$age[which (d_shellppl$anonymeID == d_shells$anonymeID[i])]
}

dat_shells <- list(
  M = nrow(d_shells),
  age = d_shells$age / mean(d_shells$age),
  tide = d_shells$tide_avg_depth
)

m_tide <- cstan(file = "models/tide_age.stan", data = dat_shells)
