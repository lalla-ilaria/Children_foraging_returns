
##################################################
#FIT MODELS AGE ONLY
##################################################
dat_shells_age <- make_list_data_age(foraging_type = "shells")
dat_traps_age <- make_list_data_age(foraging_type = "traps")

m_shell_age <- cstan( file= "models/1_shell_age.stan" , data=dat_shells_age , chains=4, cores = 4, iter = 2000 )

m_trap_age <- cstan( file= "models/1_trap_age_poisson.stan" , data=dat_traps_age , chains=4, cores = 4, iter = 2000 )

prec_shells<-precis(m_shell_age,dept=3,prob=0.95)
prec_traps<-precis(m_trap_age,dept=3,prob=0.95)
write.csv(prec_shells,"4_outcomes/model_fit/shells_prec_age.csv")
write.csv(prec_traps,"4_outcomes/model_fit/traps_prec_age.csv")
post_s<-extract.samples(m_shell_age)
post_t<-extract.samples(m_trap_age)
save(post_s, file = "4_outcomes/model_fit/post_s_age.rda")
save(post_t, file = "4_outcomes/model_fit/post_t_age.rda")



##################################################
#FIT MODELS all data
##################################################
#load data
dat_shells_all <- make_list_data_all(foraging_type = "shells")
dat_traps_all <- make_list_data_all(foraging_type = "traps")

#remove NA values - stan does not accept them - subsitute with ridiculous value
dat_shells_all[["answers_f"]][is.na(dat_shells_all[["answers_f"]])] <- -999
dat_shells_all[["answers_q"]][is.na(dat_shells_all[["answers_q"]])] <- -999
dat_shells_all[["answers_r"]][is.na(dat_shells_all[["answers_r"]])] <- -999
dat_shells_all[["knowledge_nit"]][is.na(dat_shells_all[["knowledge_nit"]])] <- -999
dat_shells_all[["height"]][is.na(dat_shells_all[["height"]])] <- -999
dat_shells_all[["grip"]][is.na(dat_shells_all[["grip"]])] <- -999


dat_traps_all[["answers_f"]][is.na(dat_traps_all[["answers_f"]])] <- -999
dat_traps_all[["answers_q"]][is.na(dat_traps_all[["answers_q"]])] <- -999
dat_traps_all[["answers_r"]][is.na(dat_traps_all[["answers_r"]])] <- -999
dat_traps_all[["knowledge_nit"]][is.na(dat_traps_all[["knowledge_nit"]])] <- -999
dat_traps_all[["height"]][is.na(dat_traps_all[["height"]])] <- -999
dat_traps_all[["grip"]][is.na(dat_traps_all[["grip"]])] <- -999

#fit models

# #run on questionnaire only
# dat_shells_all[["answers"]] <- dat_shells_all[["answers_q"]]
# dat_traps_all[["answers"]] <- dat_traps_all[["answers_q"]]
# dat_shells_all[["Q"]] <- dat_shells_all[["Q_n"]]
# dat_traps_all[["Q"]] <- dat_traps_all[["Q_n"]]
# 
# m_shells_all <- cstan( file= "models/2_shells_all.stan" , data=dat_shells_all , 
#                        chains=3, cores = 3, iter = 2000 )
# 
# m_traps_all <- cstan( file= "models/2_traps_all_poisson.stan" , data=dat_traps_all , 
#                       chains=3, cores = 3, iter = 2000 )
# #save  outputs
# prec_shellsq<-precis(m_shells_all,dept=3,prob=0.95)
# prec_trapsq<-precis(m_traps_all,dept=3,prob=0.95)
# write.csv(prec_shells,"4_outcomes/model_fit/shells_prec_onlyq.csv")
# write.csv(prec_traps,"4_outcomes/model_fit/traps_prec_onlyq.csv")
# post_sq<-extract.samples(m_shells_all)
# post_tq<-extract.samples(m_traps_all)
# save(post_s, file = "4_outcomes/model_fit/post_s_all_onlyq.rda")
# save(post_t, file = "4_outcomes/model_fit/post_t_all_onlyq.rda")

#run with all types of knowledge combined
m_shells_allk <- cstan( file= "models/2_shells_allk.stan" , data=dat_shells_all , 
                       chains=4, cores = 4, iter = 2000 )

m_traps_allk <- cstan( file= "models/2_traps_allk_poisson.stan" , data=dat_traps_all , 
                      chains=4, cores = 4, iter = 2000 )

#save  outputs
prec_shellsk<-precis(m_shells_allk,dept=3,prob=0.95)
prec_trapsk<-precis(m_traps_allk,dept=3,prob=0.95)
write.csv(prec_shells,"4_outcomes/model_fit/shells_prec_allk.csv")
write.csv(prec_traps,"4_outcomes/model_fit/traps_prec_allk.csv")
post_s_allk<-extract.samples(m_shells_allk)
post_t_allk<-extract.samples(m_traps_allk)
save(post_s_allk, file = "4_outcomes/model_fit/post_s_allk.rda")
save(post_t_allk, file = "4_outcomes/model_fit/post_t_allk.rda")

#run with traps as bernoulli only
dat_traps_age <- make_list_data_age(foraging_type = "traps")
dat_traps_age[["success"]] <- ifelse(dat_traps_age[["success"]] >= 1, 1, 0)
m_traps_age_bern <- cstan( file= "models/1_trap_age_bernoulli.stan" , data=dat_traps_age , 
                       chains=4, cores = 4, iter = 2000 )
#save  outputs
prec_traps_bern<-precis(m_traps_age_bern,dept=3,prob=0.95)
write.csv(prec_traps,"4_outcomes/model_fit/traps_prec_bern.csv")
post_t_bern<-extract.samples(m_traps_age_bern)
save(post_t_bern, file = "4_outcomes/model_fit/post_t_bern.rda")

#run complete cases only
#load data
dat_shells_cc <- make_list_data_complete_cases(foraging_type = "shells")
dat_traps_cc <- make_list_data_complete_cases(foraging_type = "traps")

m_shells_all_cc <- cstan( file= "models/2_shells_allk.stan" , data=dat_shells_cc , 
                        chains=4, cores = 4, iter = 2000 )

m_traps_all_cc <- cstan( file= "models/2_traps_allk_poisson.stan" , data=dat_traps_cc , 
                       chains=4, cores = 4, iter = 2000 )

#save  outputs
prec_shells_cc<-precis(m_shells_all_cc,dept=3,prob=0.95)
prec_traps_cc<-precis(m_traps_all_cc,dept=3,prob=0.95)
write.csv(prec_shells_cc,"4_outcomes/model_fit/shells_prec_cc.csv")
write.csv(prec_traps_cc,"4_outcomes/model_fit/traps_prec_cc.csv")
post_s_all_cc<-extract.samples(m_shells_all_cc)
post_t_all_cc<-extract.samples(m_traps_all_cc)
save(post_s_all_cc, file = "4_outcomes/model_fit/post_s_all_cc.rda")
save(post_t_all_cc, file = "4_outcomes/model_fit/post_t_all_cc.rda")




