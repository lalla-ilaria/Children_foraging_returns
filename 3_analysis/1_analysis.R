# library(rethinking)
# library(dplyr)
# library(rlist)
# real_data <- list.load("2_data_preparation/processed_data.RData")

##################################################
#FIT MODELS AGE ONLY
##################################################
dat_shells_age <- make_list_data_age(foraging_type = "shells")
dat_traps_age <- make_list_data_age(foraging_type = "traps")

m_shell_age <- cstan( file= "models/1_shell_age.stan" , data=dat_shells_age , chains=3, cores = 3, iter = 2000 )

m_trap_age <- cstan( file= "models/1_trap_age.stan" , data=dat_traps_age , chains=3, cores = 3, iter = 2000 )

prec_shells<-precis(m_shells_age,dept=3,prob=0.95)
prec_traps<-precis(m_traps_age,dept=3,prob=0.95)
write.csv(prec_shells,"4_outcomes/model_fit/shells_prec_age.csv")
write.csv(prec_traps,"4_outcomes/model_fit/traps_prec_age.csv")
post_s<-extract.samples(m_shells_age)
post_t<-extract.samples(m_traps_age)
save(post_s, file = "4_outcomes/model_fit/post_s_age.rda")
save(post_t, file = "4_outcomes/model_fit/post_t_age.rda")



##################################################
#FIT MODELS all data
##################################################
#load data
dat_shells_all <- make_list_data_all(foraging_type = "shells")
dat_traps_all <- make_list_data_all(foraging_type = "traps")

#remove NA values - stan does not accept them - subsitute with ridiculous value
dat_shells_all[["answers"]][is.na(dat_shells_all[["answers"]])] <- -999
#dat_shells_all[["knowledge"]][is.na(dat_shells_all[["knowledge"]])] <- -999
dat_shells_all[["height"]][is.na(dat_shells_all[["height"]])] <- -999
dat_shells_all[["grip"]][is.na(dat_shells_all[["grip"]])] <- -999


dat_traps_all[["answers"]][is.na(dat_traps_all[["answers"]])] <- -999
#dat_traps_all[["knowledge"]][is.na(dat_traps_all[["knowledge"]])] <- -999
dat_traps_all[["height"]][is.na(dat_traps_all[["height"]])] <- -999
dat_traps_all[["grip"]][is.na(dat_traps_all[["grip"]])] <- -999

#fit models
m_shells_all <- cstan( file= "models/9_shells_estimate_all_negative_effects.stan" , data=dat_shells_all , 
                       chains=3, cores = 3, iter = 2000 )

m_traps_all <- cstan( file= "models/9_traps_estimate_all_negative_effects.stan" , data=dat_traps , 
                      chains=3, cores = 3, iter = 2000 )

#save  outputs
prec_shells<-precis(m_shells_all,dept=3,prob=0.95)
prec_traps<-precis(m_traps_all,dept=3,prob=0.95)
write.csv(prec_shells,"4_outcomes/model_fit/shells_prec.csv")
write.csv(prec_traps,"4_outcomes/model_fit/traps_prec.csv")
post_s<-extract.samples(m_shells_all)
post_t<-extract.samples(m_traps_all)
save(post_s, file = "4_outcomes/model_fit/post_s.rda")
save(post_t, file = "4_outcomes/model_fit/post_t.rda")




##################
#TIDE MODEL
#################
dat_tides <- dat_shells_age [ c("M", "ID_i", "tide", "age")]

dat_tides$age <- d_tides$age[d_tides$ID_i]

m_tide <- cstan(file = "models/tide_age.stan", data = dat_tides)








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
