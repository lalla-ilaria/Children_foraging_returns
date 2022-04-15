#Original plots for the poster have been done with models different from those currently in the folder
#the models have the following structure
               #  \phi = \iota_i * (1-exp(-\beta_1 * age_i  )) ^ {\gamma_1} * \\
               # (1-exp(-\beta_2 * knowledge_i  )) ^ {\gamma_2} * \\
               # (1-exp(-\beta_3 * body_i  )) ^ {\gamma_3} \\
#current models use a different parameterization but produce a qualitatively similar result
#code for obtaining parallel result from the new models can be found in "1_simulation/mix cobb douglas.R"
library(rethinking)
library(rlist)
source("1_simulation/1_simulation.R")
real_data <- list.load("2_data_preparation/processed_data.RData")
trapcol <- "#5ca81e"
shellcol <-  "#ea5914"
seq_trait <- seq(0,3,0.01)

##########################################################################
#PREPARE REAL DATA
##########################################################################
#Real data
dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$knowledge),]
dc_shellppl <- dc_shellppl[complete.cases(dc_shellppl$height),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]
dat_shells <- list(
  N = nrow(dc_shellppl),
  M = nrow(dc_shells),
  A = dc_shellppl$age[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$age),
  K = dc_shellppl$knowledge[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$knowledge),
  B = dc_shellppl$height[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$height),
  R = as.numeric(dc_shells$returns)/1000,
  L = dc_shells$lenght_min/mean(dc_shells$lenght_min),
  ID_ind= as.integer(as.factor(as.character(dc_shells$anonymeID)))
)

dc_trapppl <- real_data$trap_ppl[complete.cases(real_data$trap_ppl$knowledge),]
dc_traps <- real_data$traps[which(real_data$traps$anonymeID %in% dc_trapppl$anonymeID),]
dat_traps <- list(
  N = nrow(dc_trapppl),
  M = nrow(dc_traps),
  A = dc_trapppl$age[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$age),
  K = dc_trapppl$knowledge[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$knowledge),
  B = dc_trapppl$height[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$height),
  S = as.numeric(dc_traps$success),
  L = dc_traps$lenght_hour/mean(dc_traps$lenght_hour),
  ID_ind= as.integer(as.factor(as.character(dc_traps$anonymeID)))
)

############################################################################
#AGE ONLY
############################################################################
#Real data

m_r <- cstan( file= "models/Returns_o.stan" , data=dat_shells , chains=3, cores = 3, control = list(adapt.delta = 0.9) )
m_s <- cstan( file= "models/Success_o.stan" , data=dat_traps , chains=3, cores = 3, control = list(adapt.delta = 0.9) )

post_r <- extract.samples(m_r)
post_s <- extract.samples(m_s)

#plot

png("../plots/age_only.png", height = 300, width = 450)
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,30), ylim = c(0,1), 
     xlab = "Age", ylab = "Proportion improvement")#, main = "Age only"
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))

for(i in 1:10){
  phi <- (1-exp(-post_r$beta_a[i] * seq_trait  )) ^ post_r$gamma_a[i]
  lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:10){
  phi <-  (1-exp(-post_s$beta_a[i] * seq_trait  )) ^ post_s$gamma_a[i]
  lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 0.3))
}
dev.off()

############################################################################
#KNOWLEDGE 
############################################################################
m_rk <- cstan( file= "models/Returns_k.stan" , data=dat_shells , chains=3, cores = 3, control = list(adapt.delta = 0.9) )
m_sk <- cstan( file= "models/Success_k.stan" , data=dat_traps , chains=3, cores = 3, control = list(adapt.delta = 0.9) )

post_r <- extract.samples(m_rk)
post_s <- extract.samples(m_sk)

#plot age
plot(NULL, xlim = c(0,30), ylim = c(0,1), 
     xlab = "Age", ylab = "Proportion improvement", main = "Age&Knowledge")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_a[i] * seq_trait  )) ^ post_r$gamma_a[i]
  lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_a[i] * seq_trait  )) ^ post_s$gamma_a[i]
  lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 0.3))
}

#plot knowledge
plot(NULL, xlim = c(0,3), ylim = c(0,1),
     xlab = "Knowledge", ylab = "Proportion improvement", main = "Age&Knowledge")
phi <- (1-exp(-median(post_r$beta_k) * seq_trait  )) ^ median(post_r$gamma_k)
lines( seq_trait ,  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_k) * seq_trait  )) ^ median(post_s$gamma_k)
lines( seq_trait ,  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_k * x )) ^ post_r$gamma_k, 0.95) )
shade(mu_phi, seq_trait, col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_k * x )) ^ post_s$gamma_k, 0.95) )
shade(mu_phi, seq_trait, col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_k * x )) ^ post_r$gamma_k, 0.5) )
shade(mu_phi, seq_trait, col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_k * x )) ^ post_s$gamma_k, 0.5) )
shade(mu_phi, seq_trait, col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_k * x )) ^ post_r$gamma_k, 0.3) )
shade(mu_phi, seq_trait, col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_k * x )) ^ post_s$gamma_k, 0.3) )
shade(mu_phi, seq_trait, col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_k[i] * seq_trait  )) ^ post_r$gamma_k[i]
  lines( seq_trait ,  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_k[i] * seq_trait  )) ^ post_s$gamma_k[i]
  lines( seq_trait ,  phi, col = col.alpha(trapcol, 0.3))
}

png("../plots/age_knowledge.png", height = 300, width = 450)
par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
#plot total age effect
plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * mean(dat_shells$K)  )) ^ median(post_r$gamma_k) 
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * mean(dat_traps$K)  )) ^ median(post_s$gamma_k) 
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * min(dat_shells$K)  )) ^ median(post_r$gamma_k)
mu_phi[2,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * max(dat_shells$K)  )) ^ median(post_r$gamma_k)
shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.2))

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * min(dat_shells$K)  )) ^ median(post_s$gamma_k)
mu_phi[2,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * max(dat_shells$K)  )) ^ median(post_s$gamma_k)
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.2))
dev.off()


############################################################################
#BODY
############################################################################
m_rb <- cstan( file= "models/Returns_b.stan" , data=dat_shells , chains=3, cores = 3, control = list(adapt.delta = 0.9) )
m_sb <- cstan( file= "models/Success_b.stan" , data=dat_traps , chains=3, cores = 3, control = list(adapt.delta = 0.9) )

post_s <- extract.samples(m_sb)
post_r <- extract.samples(m_rb)

#plot age
plot(NULL, xlim = c(0,30), ylim = c(0,1), 
     xlab = "Age", ylab = "Proportion improvement", main = "Age&Height")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_a[i] * seq_trait  )) ^ post_r$gamma_a[i] 
  lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_a[i] * seq_trait  )) ^ post_s$gamma_a[i]
  lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 0.3))
}

#plot body
plot(NULL, xlim = c(0,160), ylim = c(0,1), 
     xlab = "Height", ylab = "Proportion improvement", main = "Age&Height")
phi <- (1-exp(-median(post_r$beta_b) * seq_trait  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$height),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_b) * seq_trait  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$height),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_b[i] * seq_trait  )) ^ post_r$gamma_b[i]
  lines( seq_trait * mean(dc_shellppl$height),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_b[i] * seq_trait  )) ^ post_s$gamma_b[i]
  lines( seq_trait * mean(dc_trapppl$height),  phi, col = col.alpha(trapcol, 0.3))
}

png("../plots/age_height.png", height = 300, width = 450)
par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
#plot total age effect
plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_b) * mean(dat_shells$B)  )) ^ median(post_r$gamma_b) 
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_b) * mean(dat_traps$B)  )) ^ median(post_s$gamma_b) 
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_b) * min(dat_shells$B)  )) ^ median(post_r$gamma_b)
mu_phi[2,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_b) * max(dat_shells$B)  )) ^ median(post_r$gamma_b)
shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.2))

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_b) * min(dat_shells$B)  )) ^ median(post_s$gamma_b)
mu_phi[2,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_b) * max(dat_shells$B)  )) ^ median(post_s$gamma_b)
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.2))
dev.off()
############################################################################
#EVERYTHING TOGETHER
############################################################################
m_ra <- cstan( file= "models/Returns_all.stan" , data=dat_shells , chains=3, cores = 3, control = list(adapt.delta = 0.9) )
m_sa <- cstan( file= "models/Success_all.stan" , data=dat_traps , chains=3, cores = 3, control = list(adapt.delta = 0.9) )
post_r <- extract.samples(m_ra)
post_s <- extract.samples(m_sa)

#plot age
plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_a[i] * seq_trait  )) ^ post_r$gamma_a[i]
  lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_a[i] * seq_trait  )) ^ post_s$gamma_a[i]
  lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 0.3))
}

#plot knowledge effect
plot(NULL, xlim = c(0,150), ylim = c(0,1), xlab = "Knowledge", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_k) * seq_trait  )) ^ median(post_r$gamma_k)
lines( seq_trait * mean(dc_shellppl$knowledge),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_k) * seq_trait  )) ^ median(post_s$gamma_k)
lines( seq_trait * mean(dc_trapppl$knowledge),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_k * x )) ^ post_r$gamma_k, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$knowledge), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_k * x )) ^ post_s$gamma_k, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$knowledge), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_k * x )) ^ post_r$gamma_k, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$knowledge), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_k * x )) ^ post_s$gamma_k, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$knowledge), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_k * x )) ^ post_r$gamma_k, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$knowledge), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_k * x )) ^ post_s$gamma_k, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$knowledge), col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_k[i] * seq_trait  )) ^ post_r$gamma_k[i]
  lines( seq_trait * mean(dc_shellppl$knowledge),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_k[i] * seq_trait  )) ^ post_s$gamma_k[i]
  lines( seq_trait * mean(dc_trapppl$knowledge),  phi, col = col.alpha(trapcol, 0.3))
}

#plot body effect
plot(NULL, xlim = c(0,160), ylim = c(0,1), xlab = "Body", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_b) * seq_trait  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$height),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_b) * seq_trait  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$height),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_b[i] * seq_trait  )) ^ post_r$gamma_b[i]
  lines( seq_trait * mean(dc_shellppl$height),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_b[i] * seq_trait  )) ^ post_s$gamma_b[i]
  lines( seq_trait * mean(dc_trapppl$height),  phi, col = col.alpha(trapcol, 0.3))
}


#plot total age effect
plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * mean(dat_shells$K)  )) ^ median(post_r$gamma_k) *
       (1-exp(-median(post_r$beta_b) * mean(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * mean(dat_traps$K)  )) ^ median(post_s$gamma_k) *
       (1-exp(-median(post_s$beta_b) * mean(dat_traps$B)  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)


plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * max(dat_shells$K)  )) ^ median(post_r$gamma_k) *
       (1-exp(-median(post_r$beta_b) * mean(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * max(dat_traps$K)  )) ^ median(post_s$gamma_k) *
       (1-exp(-median(post_s$beta_b) * mean(dat_traps$B)  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * min(dat_shells$K)  )) ^ median(post_r$gamma_k) *
       (1-exp(-median(post_r$beta_b) * mean(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * min(dat_traps$K)  )) ^ median(post_s$gamma_k) *
       (1-exp(-median(post_s$beta_b) * mean(dat_traps$B)  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)


plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * mean(dat_shells$K)  )) ^ median(post_r$gamma_k) *
       (1-exp(-median(post_r$beta_b) * max(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * mean(dat_traps$K)  )) ^ median(post_s$gamma_k) *
       (1-exp(-median(post_s$beta_b) * max(dat_traps$B)  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * mean(dat_shells$K)  )) ^ median(post_r$gamma_k) *
       (1-exp(-median(post_r$beta_b) * min(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       (1-exp(-median(post_s$beta_k) * mean(dat_traps$K)  )) ^ median(post_s$gamma_k) *
       (1-exp(-median(post_s$beta_b) * min(dat_traps$B)  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.1))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.95) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.1))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.5) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.15))

mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_b * x )) ^ post_r$gamma_b, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(shellcol, 0.15))
mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_b * x )) ^ post_s$gamma_b, 0.3) )
shade(mu_phi, seq_trait* mean(dc_trapppl$height), col = col.alpha(trapcol, 0.15))


for(i in 1:30){
  phi <- (1-exp(-post_r$beta_b[i] * seq_trait  )) ^ post_r$gamma_b[i]
  lines( seq_trait * mean(dc_shellppl$height),  phi, col = col.alpha(shellcol, 0.3))
}
for(i in 1:30){
  phi <-  (1-exp(-post_s$beta_b[i] * seq_trait  )) ^ post_s$gamma_b[i]
  lines( seq_trait * mean(dc_trapppl$height),  phi, col = col.alpha(trapcol, 0.3))
}


#plot compared min knowledge and body vs max effect - think of what it means - I wanted things to add up but they obv don't right?
plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * min(dat_shells$K)  )) ^ median(post_r$gamma_k)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha("orange", 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
        (1-exp(-median(post_s$beta_k) * min(dat_traps$K)  )) ^ median(post_s$gamma_k)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha("lightgreen", 1), lwd = 2)

phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_b) * min(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha("firebrick1", 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
        (1-exp(-median(post_s$beta_b) * min(dat_traps$B)  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha("darkgreen", 1), lwd = 2)

phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_k) * max(dat_shells$K)  )) ^ median(post_r$gamma_k)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha("orange", 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
        (1-exp(-median(post_s$beta_k) * max(dat_traps$K)  )) ^ median(post_s$gamma_k)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha("lightgreen", 1), lwd = 2)

phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (1-exp(-median(post_r$beta_b) * max(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha("firebrick1", 1), lwd = 2)
phi <-  (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
        (1-exp(-median(post_s$beta_b) * max(dat_traps$B)  )) ^ median(post_s$gamma_b)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha("darkgreen", 1), lwd = 2)

phi <- median(post_r$tau_a)*(1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       median(post_r$tau_b)*(1-exp(-median(post_r$beta_b) * max(dat_shells$B)  )) ^ median(post_r$gamma_b)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha("firebrick1", 1), lwd = 2)

plot(NULL, xlim = c(0,3), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a)
lines( seq_trait ,  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_r$beta_k) * seq_trait  )) ^ median(post_r$gamma_k)
lines( seq_trait ,  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_r$beta_b) * seq_trait  )) ^ median(post_r$gamma_b)
lines( seq_trait ,  phi, col = col.alpha(shellcol, 1), lwd = 2)
