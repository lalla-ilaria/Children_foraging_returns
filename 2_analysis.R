
rm(list = ls())

source("project_support.R")

processed_data <- list.load("data/processed_data.RData")

##################################################
#FIT MODELS AGE ONLY
##################################################
dat_shells_age <- make_list_data_age(data = processed_data, foraging_type = "shells")
dat_traps_age <- make_list_data_age(data = processed_data, foraging_type = "traps")

m_shells_age <- cstan( file= "models/1_shell_age.stan" , data=dat_shells_age , chains=3, cores = 3, iter = 2000 )

m_traps_age <- cstan( file= "models/1_trap_age.stan" , data=dat_traps_age , chains=3, cores = 3, iter = 2000 )

prec_shells<-precis(m_shells_age,dept=3,prob=0.95)
prec_traps<-precis(m_traps_age,dept=3,prob=0.95)
write.csv(prec_shells,"model_fits/shells_prec_age.csv")
write.csv(prec_traps,"model_fits/traps_prec_age.csv")
post_s<-extract.samples(m_shells_age)
post_t<-extract.samples(m_traps_age)
save(post_s, file = "model_fits/post_s_age.rda")
save(post_t, file = "model_fits/post_t_age.rda")



##################################################
#FIT MODELS all data
##################################################
#load data
dat_shells_all <- make_list_data_all(data = processed_data, foraging_type = "shells")
dat_traps_all <- make_list_data_all(data = processed_data, foraging_type = "traps")

#remove NA values - stan does not accept them - subsitute with ridiculous value
dat_shells_all[["answers"]][is.na(dat_shells_all[["answers"]])] <- -999
dat_shells_all[["knowledge_nit"]][is.na(dat_shells_all[["knowledge_nit"]])] <- -999
dat_shells_all[["height"]][is.na(dat_shells_all[["height"]])] <- -999
dat_shells_all[["grip"]][is.na(dat_shells_all[["grip"]])] <- -999


dat_traps_all[["answers"]][is.na(dat_traps_all[["answers"]])] <- -999
dat_traps_all[["knowledge_nit"]][is.na(dat_traps_all[["knowledge_nit"]])] <- -999
dat_traps_all[["height"]][is.na(dat_traps_all[["height"]])] <- -999
dat_traps_all[["grip"]][is.na(dat_traps_all[["grip"]])] <- -999

#fit models
m_shells_all <- cstan( file= "models/9_shells_estimate_all_negative_effects.stan" , data=dat_shells_all , 
                       chains=3, cores = 3, iter = 2000 ) # ~90 minutes

m_traps_all <- cstan( file= "models/9_traps_estimate_all_negative_effects.stan" , data=dat_traps_all, 
                      chains=3, cores = 3, iter = 2000 ) #  # ~100 minutes

#save  outputs
prec_shells<-precis(m_shells_all,dept=3,prob=0.95)
prec_traps<-precis(m_traps_all,dept=3,prob=0.95)
write.csv(prec_shells,"model_fits/shells_prec.csv")
write.csv(prec_traps,"model_fits/traps_prec.csv")
post_s<-extract.samples(m_shells_all)
post_t<-extract.samples(m_traps_all)
save(post_s, file = "model_fits/post_s.rda")
save(post_t, file = "model_fits/post_t.rda")
