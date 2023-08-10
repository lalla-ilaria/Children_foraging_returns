
#load already generated quantities
generated_quantities <- list.load("2_data_preparation/generated_quantities.RData")


########################
#Proportion missing data
########################

shell_ppl <- real_data$shell_ppl[real_data$shell_ppl$anonymeID %in% real_data$shells$anonymeID,]
trap_ppl <- real_data$trap_ppl[real_data$trap_ppl$anonymeID %in% real_data$traps$anonymeID,]

nrow(shell_ppl)
sum(is.na(shell_ppl$height))
sum(is.na(shell_ppl$grip))
sum(is.na(shell_ppl$knowledge))

nrow(trap_ppl)
sum(is.na(trap_ppl$height))
sum(is.na(trap_ppl$grip))
sum(is.na(trap_ppl$knowledge))


#########################
#AGE ONLY
#########################
dat_shells <- make_list_data_age(foraging_type = "shells")
dat_traps <- make_list_data_age(foraging_type = "traps")

post_s <- extract.samples(m_shell_age)
post_t <- extract.samples(m_trap_age)

#generated data

#proportion max foraging by age 10
p_max_for_10 <- matrix(NA, 
                       nrow = 3, ncol = 2, 
                       dimnames = list(c("median", "5%PI", "94%PI"), 
                                       c("shells","traps")))
p_max_for_10[1,1] <- median((1-exp(-post_s$beta * 10/mean_age_shells  )) ^ median(post_s$gamma))
p_max_for_10[2:3,1] <- PI(((1-exp(-post_s$beta * 10/mean_age_shells  )) ^ post_s$gamma))

p_max_for_10[1,2] <- (1-exp(-median(post_t$beta) * 10/mean_age_traps  )) ^ median(post_t$gamma)
p_max_for_10[2:3,2] <- PI((1-exp(-post_t$beta * 10/mean_age_traps  )) ^ PI(post_t$gamma))
#WHY MEDIAN IS NOT BETWEEN PI????




######################################
# TIDE LEVELS
######################################
dat_shells <- make_list_data_all(foraging_type = "shells")

dat_tides <- dat_shells [ c("M", "ID_i", "tide", "age", "sex")]
dat_tides$age <- d_tides$age[d_tides$ID_i]
dat_tides$sex <- d_tides$sex[d_tides$ID_i]
dat_tides$age <- dat_tides$age * mean_age_shells

dat_tides <- as.data.frame(dat_tides)
dat_tides <- dat_tides [-which(dat_tides$age >=15 & dat_tides$sex == "m"),]
dat_tides %>% filter(age >=20 ) %>% summarize( Mean = mean(tide))
dat_tides %>% filter(age <=19 ) %>% summarize( Mean = mean(tide))

dat_tides %>% filter(tide >0 ) %>% summarize( Mean = mean(age))
dat_tides %>% filter(tide <=0 ) %>% summarize( Mean = mean(age))

###############
#shells calculation
dat_shells <- make_list_data_all(foraging_type = "shells")
post_s <- extract.samples(m_shells_all)

phi <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* min(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(180/mean(real_data$shells$lenght_min))) +
  post_s$tau* 0
kg_shells <- exp (log(post_s$alpha) + phi_min + psi +
                    (post_s$sigma^2 /2))
mean(kg_shells)
PI(kg_shells)

##########
#SAVE DATA
##########
generated_data <- list(p_max_for_10 = p_max_for_10)

list.save(generated_data, "4_outcomes/generated_data.Rdata")
