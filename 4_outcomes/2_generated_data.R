
#load already generated quantities
generated_quantities <- list.load("2_data_preparation/generated_quantities.RData")

#########################
#AGE ONLY
#########################
dat_shells <- make_list_data_age(foraging_type = "shells")
dat_traps <- make_list_data_age(foraging_type = "traps")

load(file = "4_outcomes/model_fit/post_s_age.rda")
load(file = "4_outcomes/model_fit/post_t_age.rda")

#generated data

#proportion max foraging by age 10
generated_quantities$shell_median_foraging_by10 <- (1-exp(-median(post_s$beta) * 10/mean_age_shells  )) ^ median(post_s$gamma)
generated_quantities$shell_5PI_foraging_by10 <-  PI((1-exp(-post_s$beta * 10/mean_age_shells  )) ^ post_s$gamma)[1]
generated_quantities$shell_94PI_foraging_by10 <-  PI((1-exp(-post_s$beta * 10/mean_age_shells  )) ^ post_s$gamma)[2]

generated_quantities$trap_median_foraging_by10 <- (1-exp(-median(post_t$beta) * 10/mean_age_traps  )) ^ median(post_t$gamma)
generated_quantities$trap_5PI_foraging_by10 <-  PI((1-exp(-post_t$beta * 10/mean_age_traps  )) ^ post_t$gamma)[1]
generated_quantities$trap_94PI_foraging_by10 <-  PI((1-exp(-post_t$beta * 10/mean_age_traps  )) ^ post_t$gamma)[2]
#WHY MEDIAN IS NOT BETWEEN PI????




######################################
# TIDE LEVELS
######################################
dat_shells <- make_list_data_all(foraging_type = "shells")

dat_tides <- dat_shells [ c("M", "ID_i", "tide", "age", "sex")]
dat_tides$age <- dat_tides$age[dat_tides$ID_i]
dat_tides$sex <- dat_tides$sex[dat_tides$ID_i]
dat_tides$age <- dat_tides$age * mean_age_shells

dat_tides <- as.data.frame(dat_tides)
generated_quantities$mean_tide_above20 <- dat_tides %>% filter(age >=20 ) %>% summarize( Mean = mean(tide)) %>% unlist(use.names = FALSE) 
generated_quantities$mean_tide_below19 <- dat_tides %>% filter(age <=19 ) %>% summarize( Mean = mean(tide)) %>% unlist(use.names = FALSE)

generated_quantities$mean_age_tide_above0 <- dat_tides %>% filter(tide >0 ) %>% summarize( Mean = mean(age)) %>% unlist(use.names = FALSE)
generated_quantities$mean_age_tide_below0 <- dat_tides %>% filter(tide <=0 ) %>% summarize( Mean = mean(age))%>% unlist(use.names = FALSE)

###############
#shells calculation
dat_shells <- make_list_data_all(foraging_type = "shells")
load(file = "4_outcomes/model_fit/post_s_all.rda")

phi <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* min(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(180/mean(real_data$shells$lenght_min))) +
  post_s$tau* 0
kg_shells <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
generated_quantities$mean_kg_shells <- mean(kg_shells)
generated_quantities$PI5_kg_shells <- PI(kg_shells)[1]
generated_quantities$PI5_kg_shells <- PI(kg_shells)[2]

##########
#SAVE DATA
##########

list.save(generated_quantities, "4_outcomes/generated_data.RData")
