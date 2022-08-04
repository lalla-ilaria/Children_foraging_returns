library(rethinking)
library(rlist)
library(tidyverse)
library(ggridges)
real_data <- list.load("2_data_preparation/processed_data.RData")
seq_trait <- seq(0,3,0.001)

#########################
#AGE ONLY
#########################
#prepare data age only
dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$age),]
dc_shell_k <- real_data$shell_k[which(rownames(real_data$shell_k) %in% dc_shellppl$anonymeID),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]

dc_shellppl$index_id <- as.integer(as.factor(dc_shellppl$anonymeID))
dc_shellppl <- dc_shellppl[order(dc_shellppl$index_id),]
dc_shells$index_id <- as.integer(as.factor(dc_shells$anonymeID))
dc_shell_k <- dc_shell_k[ order(as.factor(row.names(dc_shell_k))), ]

dc_trapppl <- real_data$trap_ppl[complete.cases(real_data$trap_ppl$age),]
dc_traps <- real_data$traps[which(real_data$traps$anonymeID %in% dc_trapppl$anonymeID),]
dc_trap_k <- real_data$trap_k[which(rownames(real_data$trap_k) %in% dc_trapppl$anonymeID),]

dc_trapppl$index_id <- as.integer(as.factor(dc_trapppl$anonymeID))
dc_trapppl <- dc_trapppl[order(dc_trapppl$index_id),]
dc_traps$index_id <- as.integer(as.factor(dc_traps$anonymeID))
dc_trap_k <- dc_trap_k[ order(as.factor(row.names(dc_trap_k))), ]

post_s <- extract.samples(m_shell_age)
post_t <- extract.samples(m_trap_age)

#generated data

#proportion max foraging by age 10
p_max_for_10 <- matrix(NA, 
                       nrow = 3, ncol = 2, 
                       dimnames = list(c("median", "5%PI", "94%PI"), 
                                       c("shells","traps")))
p_max_for_10[1,1] <- median((1-exp(-post_s$beta * 10/mean(dc_shellppl$age)  )) ^ median(post_s$gamma))
p_max_for_10[2:3,1] <- PI(((1-exp(-post_s$beta * 10/mean(dc_shellppl$age)  )) ^ post_s$gamma))

p_max_for_10[1,2] <- (1-exp(-median(post_t$beta) * 10/mean(dc_trapppl$age)  )) ^ median(post_t$gamma)
p_max_for_10[2:3,2] <- PI((1-exp(-post_t$beta * 10/mean(dc_trapppl$age)  )) ^ PI(post_t$gamma))
#WHY MEDIAN IS NOT BETWEEN PI????

generated_data <- list(p_max_for_10 = p_max_for_10)

list.save(generated_data, "4_outcomes/generated_data.Rdata")
