library(rethinking)
library(rlist)
library(tidyverse)
library(ggridges)
real_data <- list.load("2_data_preparation/processed_data.RData")
trapcol <- "#52b502" #"#5ca81e"
shellcol <-  "#ea5914"
othercol <- "grey30"
seq_trait <- seq(0,3,0.001)
age_plot <- 40

#subset data
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

##############################
#fit to data - age only
##############################
#plot fit to data with age only
post_s <- extract.samples(m_shell_age)
post_t <- extract.samples(m_trap_age)

png("../plots/fit_to_data.png", height = 10, width = 16, units = "cm", res = 500)
  par(mfrow = c(1,2),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
  #shells
  phi <-  exp(mean(post_s$iota) ) * (
      (1-exp(- mean(post_s$beta) * seq_trait  )) ^ mean(post_s$gamma))
  psi <-   (mean(dat_shells$L)) ^ mean(post_s$xi) * 
      exp(mean(dat_shells$tide) * mean(post_s$tau))
  R <- exp (  log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))
  samp_data <- rlnorm(length(seq_trait),  
                      log(mean(post_s$alpha) * phi * psi), 
                      mean(post_s$sigma))
  
  plot(jitter(seq_trait) * mean(dc_shellppl$age), samp_data, 
       xlim = c(0,age_plot), ylim = c(0, max(dat_shells$R)+1), 
       xlab = "Age", ylab = "kg shellfish",
       pch = 16, col = col.alpha("orange", 0.2))
  for(i in 1:150){
    phi <-  exp(apply(post_s$iota,1,mean )[i] ) * (
      (1-exp(- post_s$beta[i] * seq_trait  )) ^ post_s$gamma[i]
    )
    psi <-   (mean(dat_shells$L)) ^ post_s$xi[i] * 
      exp(mean(dat_shells$tide) * post_s$tau[i])
    R <- exp (  log(post_s$alpha[i] * phi * psi) + 
                  (post_s$sigma[i]^2 /2))
    lines( seq_trait * mean(dc_shellppl$age),  R, 
           col = col.alpha(shellcol, 0.2), lwd = 1)
  }
  points(jitter(dc_shellppl$age[dc_shells$index_id], amount = 0.4), 
         as.numeric(dc_shells$returns)/1000,
         pch = 16, cex = 0.8, col = col.alpha(othercol, 0.4))

  #traps 
  phi <-  exp(mean(post_t$iota) ) * (
      (1-exp(- mean(post_t$beta) * seq_trait  )) ^ mean(post_t$gamma))
  psi <-   (mean(dat_traps$L)) ^ mean(post_t$xi)
  p <- 1 - exp ( - mean(post_t$alpha) * phi * psi)
  samp_data <- rbern(length(seq_trait),  p)
  
  #NB making plot with 13212 only (one of best hunters) to make plot with optimal situation to raise curve off zero
  plot(jitter(seq_trait) * mean(dc_trapppl$age), samp_data, 
       xlab = "Age", ylab = "p capture",
       xlim = c(0,age_plot), ylim = c(0, 1), 
       pch = 16, col = col.alpha("lawngreen", 0.3))
  for(i in 1:150){
    phi <-  exp(post_t$iota[i,4] ) * (
      (1-exp(- post_t$beta[i] * seq_trait  )) ^ post_t$gamma[i]
    )
    psi <-   (min(dat_traps$L)) ^ post_t$xi[i]
    p <- 1 - exp ( - post_t$alpha[i] * phi * psi)
    lines( seq_trait * mean(dc_trapppl$age),  p, 
           col = col.alpha(trapcol, 0.2), lwd = 1)
  }
  points(jitter(dc_trapppl$age[dc_traps$index_id], amount = 0.5), 
         jitter(dc_traps$success, amount = 0.02), 
    pch = 16, cex = ifelse(dc_traps$success == 1, 0.8, 0.6), 
    col = col.alpha(othercol, ifelse(dc_traps$success == 1, 0.4, 0.1)))
  #donow6 how to show better? jitter y expand axis/try scale points
dev.off()
#######################################
#age variation
#######################################
png("../plots/phi_age_only.png", height = 10, width = 16, units = "cm", res = 500)
  par(mfrow = c(1,2),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
  plot(NULL, xlim = c(0,age_plot), ylim = c(0,1), 
       xlab = "Age", ylab = "proportion max foraging")#, main = "Age only"
  phi <- (1-exp(-median(post_s$beta) * seq_trait  )) ^ median(post_s$gamma)
  lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta * x )) ^ post_s$gamma, 0.95) )
  shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.1))
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta * x )) ^ post_s$gamma, 0.5) )
  shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta * x )) ^ post_s$gamma, 0.3) )
  shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
  for(i in 1:30){
    phi <- (1-exp(-post_s$beta[i] * seq_trait  )) ^ post_s$gamma[i]
    lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 0.3))
  }
  
  plot(NULL, xlim = c(0,age_plot), ylim = c(0,1), 
       xlab = "Age", ylab = "proportion max foraging")#, main = "Age only"
  phi <-  (1-exp(-median(post_t$beta) * seq_trait  )) ^ median(post_t$gamma)
  lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_t$beta * x )) ^ post_t$gamma, 0.95) )
  shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_t$beta * x )) ^ post_t$gamma, 0.5) )
  shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_t$beta * x )) ^ post_t$gamma, 0.3) )
  shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))
  for(i in 1:30){
    phi <-  (1-exp(-post_t$beta[i] * seq_trait  )) ^ post_t$gamma[i]
    lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 0.3))
  }
dev.off()
  
#####################################
#effect of knowledge and other stuff
#####################################
#subset data
dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$knowledge),]
dc_shellppl <- dc_shellppl[complete.cases(dc_shellppl$height),]
dc_shell_k <- real_data$shell_k[which(rownames(real_data$shell_k) %in% dc_shellppl$anonymeID),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]

dc_trapppl <- real_data$trap_ppl[complete.cases(real_data$trap_ppl$knowledge),]
dc_traps <- real_data$traps[which(real_data$traps$anonymeID %in% dc_trapppl$anonymeID),]
dc_trap_k <- real_data$trap_k[which(rownames(real_data$trap_k) %in% dc_trapppl$anonymeID),]

#make indexes
dc_shellppl$index_id <- as.integer(as.factor(dc_shellppl$anonymeID))
dc_shellppl <- dc_shellppl[order(dc_shellppl$index_id),]
dc_shells$index_id <- as.integer(as.factor(dc_shells$anonymeID))
dc_shell_k <- dc_shell_k[ order(as.factor(row.names(dc_shell_k))), ]

dc_trapppl$index_id <- as.integer(as.factor(dc_trapppl$anonymeID))
dc_trapppl <- dc_trapppl[order(dc_trapppl$index_id),]
dc_traps$index_id <- as.integer(as.factor(dc_traps$anonymeID))
dc_trap_k <- dc_trap_k[ order(as.factor(row.names(dc_trap_k))), ]

#extract samples
post_s <- extract.samples(m_shell_all)
post_t <- extract.samples(m_trap_all)
#make this code legible :o
#shells
#knowledge
phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      min(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
min_k <- exp (log(mean(post_s$alpha) * phi * psi) + 
              (mean(post_s$sigma)^2 /2))
#rlnorm(length(seq_trait),  #change to expectation, no simulation
#                      log(post_s$alpha * phi * psi), 
#                      post_s$sigma)

phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      max(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
max_k <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))
diff_k <- max_k - min_k
plot(density(diff_k))

#grip
phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      min(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
min_g <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))

phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      max(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
max_g <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))
diff_g <- max_g - min_g
plot(density(diff_g))

#height
phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      min(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
min_h <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))

phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      max(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
max_h <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))
diff_h <- max_h - min_h
plot(density(diff_h))

#tide
phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(min(dat_shells$tide) * post_s$tau)
min_t <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))

phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (mean(dat_shells$L)) ^ post_s$xi * 
    exp(max(dat_shells$tide) * post_s$tau)
max_t <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))
diff_t <- max_t - min_t
plot(density(diff_t))

#duration
phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (min(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
min_l <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))

phi <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
psi <-   (max(dat_shells$L)) ^ post_s$xi * 
    exp(mean(dat_shells$tide) * post_s$tau)
max_l <- exp (log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))
diff_l <- max_l - min_l
plot(density(diff_l))


diffs <- data.frame(variable = c( rep("knowledge", length(diff_k)),
                               rep("grip", length(diff_g)),
                               rep("height", length(diff_h)),
                               rep("duration", length(diff_l)),
                               rep("tide", length(diff_t))),
                    kg_shellfish =  c(diff_k, diff_g, diff_h, diff_l, -diff_t)
                    )

diffs <- diffs %>% 
  mutate( variable = fct_reorder(.f = variable, .x = kg_shellfish, .fun = mean))
ggplot(diffs, aes(x = kg_shellfish, y = variable)) + 
  geom_vline(xintercept = 0, colour = col.alpha("grey40", 0.3))+
  geom_density_ridges(fill = col.alpha("cornflowerblue", 0.3), 
                      colour = "cornflowerblue") +
  labs(x = "max-min kg shellfish", y = "")+
  xlim(-0.5, 8)+
  theme_classic()


#traps
#knowledge
phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      min(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (mean(dat_traps$L)) ^ post_t$xi 
min_k <- 1 - exp ( - post_t$alpha * phi * psi)

  
phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      max(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (mean(dat_traps$L)) ^ post_t$xi 
max_k <- 1 - exp ( - post_t$alpha * phi * psi)

diff_k <- max_k - min_k

#grip
phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      min(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (mean(dat_traps$L)) ^ post_t$xi 
min_g <- 1 - exp ( - post_t$alpha * phi * psi)

phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      max(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (mean(dat_traps$L)) ^ post_t$xi 
max_g <- 1 - exp ( - post_t$alpha * phi * psi)
diff_g <- max_g - min_g

#height
phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      min(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (mean(dat_traps$L)) ^ post_t$xi 
min_h <- 1 - exp ( - post_t$alpha * phi * psi)

phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      max(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (mean(dat_traps$L)) ^ post_t$xi
max_h <- 1 - exp ( - post_t$alpha * phi * psi)
diff_h <- max_h - min_h

#duration
phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (min(dat_traps$L)) ^ post_t$xi 
min_l <- 1 - exp ( - post_t$alpha * phi * psi)

phi <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
psi <-   (max(dat_traps$L)) ^ post_t$xi 
max_l <- 1 - exp ( - post_t$alpha * phi * psi)
diff_l <- max_l - min_l

diffs <- data.frame(variable = c( rep("knowledge", length(diff_k)),
                               rep("grip", length(diff_g)),
                               rep("height", length(diff_h)),
                               rep("duration", length(diff_l))),
                    p_capture =  c(diff_k, diff_g, diff_h, diff_l)
                    )

diffs <- diffs %>% 
  mutate( variable = fct_reorder(.f = variable, .x = p_capture, .fun = mean))
ggplot(diffs, aes(x = p_capture, y = variable)) + 
  geom_vline(xintercept = 0, colour = col.alpha("grey40", 0.3))+
  geom_density_ridges(fill = col.alpha("cornflowerblue", 0.3), 
                      colour = "cornflowerblue") +
  labs(x = "max-min p of trap capturing", y = "")+
  theme_classic()

#######################################
#effect of variables in shells vs traps
#######################################
#shells
#knowledge
phi_min_k <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      min(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )

phi_max_k <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      max(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
phi_k <- phi_max_k - phi_min_k

#grip
phi_min_g <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      min(dat_shells$G) ^ post_s$theta_g 
  )

phi_max_g <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      mean(dat_shells$H) ^ post_s$eta_h *
      max(dat_shells$G) ^ post_s$theta_g 
  )
phi_g <- phi_max_g - phi_min_g

#height
phi_min_h <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      min(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )

phi_max_h <-  exp(apply(post_s$iota,1,mean ) ) * (
    (1-exp(- post_s$beta * 20  )) ^ post_s$gamma *
      mean(apply(post_s$K, 2, mean)) ^ post_s$zeta_k  *
      max(dat_shells$H) ^ post_s$eta_h *
      mean(dat_shells$G) ^ post_s$theta_g 
  )
phi_h <- phi_max_h - phi_min_h
diffs_shells <- data.frame(variable = c( rep("knowledge", length(phi_k)),
                               rep("grip", length(phi_g)),
                               rep("height", length(phi_h))),
                    phi =  c(phi_k, phi_g, phi_h),
                    type = "shells"
                    )


#traps
#knowledge
phi_min_k <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      min(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
  
phi_max_k <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      max(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
phi_k <- phi_max_k - phi_min_k

#grip
phi_min_g <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      min(dat_traps$G) ^ post_t$theta_g 
  )

phi_max_g <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      mean(dat_traps$H) ^ post_t$eta_h *
      max(dat_traps$G) ^ post_t$theta_g 
  )
phi_g <- phi_max_g - phi_min_g

#height
phi_min_h <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      min(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )

phi_max_h <-  exp(apply(post_t$iota,1,mean ) ) * (
    (1-exp(- post_t$beta * 20  )) ^ post_t$gamma *
      mean(apply(post_t$K, 2, mean)) ^ post_t$zeta_k  *
      max(dat_traps$H) ^ post_t$eta_h *
      mean(dat_traps$G) ^ post_t$theta_g 
  )
phi_h <- phi_max_h - phi_min_h

diffs_traps <- data.frame(variable = c( rep("knowledge_t", length(phi_k)),
                               rep("grip_t", length(phi_g)),
                               rep("height_t", length(phi_h))),
                    phi =  c(phi_k, phi_g, phi_h),
                    type = "traps"
                    )
diffs <- rbind(diffs_shells, diffs_traps)
# diffs <- diffs %>% 
#   mutate( variable = fct_reorder(.f = variable, .x = p, .fun = mean))
ggplot(diffs, aes(x = phi, y = variable, fill = type)) + 
  geom_vline(xintercept = 0, colour = col.alpha("grey40", 0.3))+
  geom_density_ridges(colour = col.alpha("cornflowerblue", 0.3)) +
  xlim(-0.1,2)+
  theme_classic()

