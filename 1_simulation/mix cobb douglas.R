library(rethinking)
library(rlist)
library(tidyverse)
source("1_simulation/1_simulation.R")
real_data <- list.load("2_data_preparation/processed_data.RData")
trapcol <- "#5ca81e"
shellcol <-  "#ea5914"
seq_trait <- seq(0,3,0.001)

#####################
nsamp <-100
#PRIOR PREDICTIVE SIMULATION
#Set trait values for age and possible times
AGE <- seq(0,3, 0.1) #trait
L <- c(0, 0.1, 1, 3) #trait
#create matrixes to save individual and trip effects
phi <- matrix(NA, nsamp, length(AGE))
psi <- rep(NA, nsamp)
#SIMULATED PRIORS
#general intercept
alpha <- rnorm(nsamp, 0, 0.1)#prior
#trip
lambda <-  rexp(nsamp,1)
#individual_age
beta <- rlnorm(nsamp, 0, 1)#prior
gamma <- rlnorm(nsamp, 0, 1)#prior
#sigma lognormal
sigma <- rexp(nsamp, 1)
#ADD KNOWLEDGE
betak <- rlnorm(nsamp, 0, 1)#prior
zeta_k <- rlnorm(nsamp, 0, 1)#prior
K <- matrix(NA, nsamp, length(AGE)) 
for(i in 1:nsamp){
  K[i,] <- betak[i]*AGE #trait
}
betab <- rlnorm(nsamp, 0, 1)#prior
eta_b <- rlnorm(nsamp, 0, 1)#prior
B <- matrix(NA, nsamp, length(AGE)) 
for(i in 1:nsamp){
  B[i,] <- betab[i]*AGE #trait
}

#plot prior predictive simulation
plot(NULL, xlim = c(0,3), ylim = c (0,3), 
     xlab = "Age", ylab = "Proportion improvement")
#calculate per sample
for(i in 1:nsamp){
  phi[i,] <- (1-exp(-beta[i] * AGE  )) ^ gamma[i] 
  lines( AGE,  phi[i,], col = col.alpha("cornflowerblue", 0.3))
}
for(i in 1:nsamp){
  phi[i,] <- (1-exp(-beta[i] * AGE  )) ^ gamma[i] *
    ( mean(K) ) ^ zeta_k[i]
  lines( AGE,  phi[i,], col = col.alpha("lightblue", 0.3))
}

lines (AGE, (1-exp(-1 * AGE  )) ^ 1 ,
       col = col.alpha("darkblue", 0.7))
lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( mean(K)  ) ^ -1,
       col = col.alpha("darkblue", 0.7))
lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( mean(K)  ) ^ 0.1,
       col = col.alpha("darkblue", 0.7))
lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( mean(K)  ) ^ 0,
       col = col.alpha("darkred", 0.7))

lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( 0.1 ) ^ 1,
       col = col.alpha("darkblue", 0.7))


#SIMULATE
#Simulate
d<-sim_data(50,300, zero = F, b_a = 1, g_a = 1, g_k = 1, g_b = 1)
d_simp <-d
par(mfrow = c(2,2),mgp = c(1.5, 0.5, 0), mar = c(2.2, 2.2, 1.5, 0.5) + 0.1)
plot(d$B, d$K)
plot(d$A[d$ID_ind], d$R)
plot(d$K[d$ID_ind], d$R)
plot(d$B[d$ID_ind], d$R)
dat <- list(
  N = d$N,
  M = d$M,
  A = d$A/mean(d$A),
  K = d$K/mean(d$K),
  B = d$B/mean(d$B),
  R = d$R,
  L = d$L/mean(d$L),
  ID_i = d$ID_ind
)


m_ro <- cstan( file= "models/Returns_o.stan" , data=dat , chains=3, cores = 3 )
post_r <- extract.samples(m_ro)
plot(precis(m_ro))

plot(d$phi, apply(post_r$phi, 2, mean))
for(i in 1:50) lines(rep(d$phi[i],2), apply(post_r$phi, 2, PI)[,i])

curve((1-exp(- post_r$beta[1] * x  )) ^ post_r$gamma[1] , xlim = c(0,3), ylim = c(0,1), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i] , col = col.alpha("grey80", 0.5), add = T)

plot(d$A[d$ID_ind], dat$R, 
     xlim = c(0,20), ylim = c(0, max(dat$R)+1), pch = 16, col = col.alpha("grey40", 0.7),
     main = "gamma = lnorm(1,1)")
for(i in 1:150){
  phi <-  exp(apply(post_r$iota,1,mean )[i] ) * (
    (1-exp(- post_r$beta[i] * seq(0,3,0.1)  )) ^ post_r$gamma[i] 
  )
  psi <-   (mean(dat$L)) ^ post_r$xi[i]
  R <-  exp(log(post_r$alpha[i] * phi * psi) + (post_r$sigma[i]^2 /2))
  samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_r$alpha[i] * phi * psi), post_r$sigma[i])
  points(jitter(seq(0,3,0.1)) * mean(d$A), samp_data, col = col.alpha("orange", 0.1), pch = 16)
  lines( seq(0,3,0.1) * mean(d$A),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}


m_ras <- cstan( file= "models/Returns_all.stan" , data=dat , chains=3, cores = 3 )
post_r <- extract.samples(m_ras)
par(mfrow = c(2,2))
plot(precis(m_ras))

plot(d$phi, apply(post_r$phi, 2, mean))
for(i in 1:50) lines(rep(d$phi[i],2), apply(post_r$phi, 2, PI)[,i])

curve((1-exp(- post_r$beta[1] * x  )) ^ post_r$gamma[1] , xlim = c(0,3), ylim = c(0,2), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i] , col = col.alpha("grey80", 0.5), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat$K) ^ post_r$zeta_k[i]  *
                         min(dat$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat$K) ^ post_r$zeta_k[i]  *
                         max(dat$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         min(dat$K) ^ post_r$zeta_k[i]  *
                         mean(dat$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         max(dat$K) ^ post_r$zeta_k[i]  *
                         mean(dat$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)
plot(d$A[d$ID_ind], dat$R, 
     xlim = c(0,20), ylim = c(0, max(dat$R)+1), pch = 16, col = col.alpha("grey40", 0.7),
     main = "alpha = norm(1,1)")
for(i in 1:150){
phi <-  exp(apply(post_r$iota,1,mean )[i] ) * (
   (1-exp(- post_r$beta[i] * seq(0,3,0.1)  )) ^ post_r$gamma[i] *
    mean(dat$K) ^ post_r$zeta_k[i]  *
    mean(dat$B) ^ post_r$eta_b[i] 
  )
psi <-   (mean(dat$L)) ^ post_r$xi[i]
R <- exp (  log(post_r$alpha[i] * phi * psi) + (post_r$sigma[i]^2 /2))
samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_r$alpha[i] * phi * psi), post_r$sigma[i])
points(jitter(seq(0,3,0.1)) * mean(d$A), samp_data, col = col.alpha("orange", 0.1), pch = 16)
lines( seq(0,3,0.1) * mean(d$A),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}


m_rool <- cstan( file= "models/Returns_all_outoflog.stan" , data=dat , chains=3, cores = 3 )
post_r <- extract.samples(m_rool)
par(mfrow = c(2,2))
plot(precis(m_rool))

plot(d$phi, apply(post_r$phi, 2, mean))
for(i in 1:50) lines(rep(d$phi[i],2), apply(post_r$phi, 2, PI)[,i])

curve(log(1-exp(- post_r$beta[1] * x  )) * post_r$gamma[1] , xlim = c(0,3), ylim = c(0,2), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i] , col = col.alpha("grey80", 0.5), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat$K) ^ post_r$zeta_k[i]  *
                         min(dat$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat$K) ^ post_r$zeta_k[i]  *
                         max(dat$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         min(dat$K) ^ post_r$zeta_k[i]  *
                         mean(dat$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         max(dat$K) ^ post_r$zeta_k[i]  *
                         mean(dat$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)
plot(d$A[d$ID_ind], dat$R, 
     xlim = c(0,20), ylim = c(0, max(dat$R)+1), pch = 16, col = col.alpha("grey40", 0.7),
     main = "alpha = lognormal(0,1)")
for(i in 1:150){
  phi <-  exp(apply(post_r$iota,1,mean )[i] ) * (
    (1-exp(- post_r$beta[i] * seq(0,3,0.1)  )) ^ post_r$gamma[i] *
      mean(dat$K) ^ post_r$zeta_k[i]  *
      mean(dat$B) ^ post_r$eta_b[i] 
  )
  psi <-   (mean(dat$L)) ^ post_r$xi[i]
  R <- exp (  post_r$alpha[i] + phi + psi + (post_r$sigma[i]^2 /2))
  samp_data <- rlnorm(length(seq(0,3,0.1)),  post_r$alpha[i] + phi + psi, post_r$sigma[i])
  points(jitter(seq(0,3,0.1)) * mean(d$A), samp_data, col = col.alpha("orange", 0.1), pch = 16)
  lines( seq(0,3,0.1) * mean(d$A),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}


d<-sim_data(50,300, zero = F, b_a = 1, g_a = 1, g_k = 1, g_b = 3, add_to_b = 3)
d_add <- d
par(mfrow = c(2,2))
plot(d$B, d$K)
plot(d$A[d$ID_ind], d$R)
plot(d$K[d$ID_ind], d$R)
plot(d$B[d$ID_ind], d$R)
dat <- list(
  N = d$N,
  M = d$M,
  A = d$A/mean(d$A),
  K = d$K/mean(d$K),
  B = d$B/mean(d$B),
  R = d$R,
  L = d$L/mean(d$L),
  ID_i = d$ID_ind
)
plot(d$B, d$K)
plot(dat$A[d$ID_ind], d$R)
plot((d$K/min(d$K))[d$ID_ind], d$R, xlim = c(0,3))
plot((d$B/min(d$B))[d$ID_ind], d$R, xlim = c(0,3))

m_rasp <- cstan( file= "models/Returns_all.stan" , data=dat , chains=3, cores = 3 )
post_r <- extract.samples(m_rasp)
plot(precis(m_rasp))

plot(d$phi, apply(post_r$phi, 2, mean))
for(i in 1:50) lines(rep(d$phi[i],2), apply(post_r$phi, 2, PI)[,i])

curve((1-exp(- post_r$beta[1] * x  )) ^ post_r$gamma[1] , xlim = c(0,3), ylim = c(0,2), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i] , col = col.alpha("grey80", 0.5), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat$K) ^ post_r$zeta_k[i]  *
                         min(dat$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat$K) ^ post_r$zeta_k[i]  *
                         max(dat$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         min(dat$K) ^ post_r$zeta_k[i]  *
                         mean(dat$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         max(dat$K) ^ post_r$zeta_k[i]  *
                         mean(dat$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)
plot(d$A[d$ID_ind], dat$R, 
     xlim = c(0,20), ylim = c(0, max(dat$R)+1), pch = 16, col = col.alpha("grey40", 0.7))
for(i in 1:150){
phi <-  exp(apply(post_r$iota,1,mean )[i] ) * (
   (1-exp(- post_r$beta[i] * seq(0,3,0.1)  )) ^ post_r$gamma[i] *
    mean(dat$K) ^ post_r$zeta_k[i]  *
    mean(dat$B) ^ post_r$eta_b[i] 
  )
psi <-   (mean(dat$L)) ^ post_r$xi[i]
R <- exp (  log(post_r$alpha[i] * phi * psi) + (post_r$sigma[i]^2 /2))
samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_r$alpha[i] * phi * psi), post_r$sigma[i])
points(jitter(seq(0,3,0.1)) * mean(d$A), samp_data, col = col.alpha("orange", 0.1), pch = 16)
lines( seq(0,3,0.1) * mean(d$A),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}

par(mfrow = c(1,1))


d<-sim_data(100,300, zero = T)
dat <- list(
  N = d$N,
  M = d$M,
  A = d$A/mean(d$A),
  K = d$K/mean(d$K),
  B = d$B/mean(d$B),
  S = d$S,
  L = d$L/mean(d$L),
  ID_ind = d$ID_ind
)
m_sas <- cstan( file= "models/Success_all.stan" , data=dat , chains=3, cores = 3 )

post_s <- extract.samples(m_sas)



plot(NULL, xlim = c(0,3), ylim = c(0,1))
for(i in 1:150){
  phi <- (1-exp(-post_s$beta[i] * seq_trait  )) ^ post_s$gamma[i]
  lines( seq_trait ,  phi, col = col.alpha("cornflowerblue", 0.4))
}


##########################################################################
#PREPARE REAL DATA
##########################################################################
source("3_analysis/1_analysis.R")


#mix cobb douglas

m_ra <- cstan( file= "models/Returns_all.stan" , data=dat_shells , chains=3, cores = 3 )
# m_rab <- cstan( file= "models/Returns_allbody.stan" , data=dat_shells , chains=3, cores = 3 )
m_sa <- cstan( file= "models/Success_all.stan" , data=dat_traps , chains=3, cores = 3 )

par(mfrow = c(2,2))
post_s <- extract.samples(m_sa)
post_r <- extract.samples(m_rairt)
plot(precis(m_ra))
plot(dat_shells$A, apply(post_r$phi, 2, mean), ylim=c(0,4))
 for(i in 1:dat_shells$N) lines(rep(dat_shells$A[i],2), apply(post_r$phi, 2, PI)[,i])
# plot(dat_shells$K, apply(post_r$phi, 2, mean), ylim=c(0,4))
# for(i in 1:dat_shells$N) lines(rep(dat_shells$K[i],2), apply(post_r$phi, 2, PI)[,i])
# plot(dat_shells$B, apply(post_r$phi, 2, mean), ylim=c(0,4))
# for(i in 1:dat_shells$N) lines(rep(dat_shells$B[i],2), apply(post_r$phi, 2, PI)[,i])

curve((1-exp(- post_r$beta[1] * x  )) ^ post_r$gamma[1] , xlim = c(0,3), ylim = c(0,2), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i] , col = col.alpha("grey80", 0.5), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat_shells$K) ^ post_r$zeta_k[i]  *
                         min(dat_shells$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         mean(dat_shells$K) ^ post_r$zeta_k[i]  *
                         max(dat_shells$B) ^ post_r$eta_b[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         min(dat_shells$K) ^ post_r$zeta_k[i]  *
                         mean(dat_shells$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i]  *
                         max(dat_shells$K) ^ post_r$zeta_k[i]  *
                         mean(dat_shells$B) ^ post_r$eta_b[i] , col = col.alpha(trapcol, 0.1), add = T)

plot(dc_shellppl$age[dat_shells$ID_i], dat_shells$R, 
     xlim = c(0,20), ylim = c(0, max(dat_shells$R)+1), pch = 16, col = col.alpha("grey40", 0.7))
for(i in 1:150){
  phi <-  exp(apply(post_r$iota,1,mean )[i] ) * (
    (1-exp(- post_r$beta[i] * seq(0,3,0.1)  )) ^ post_r$gamma[i] *
      mean(dat_shells$K) ^ post_r$zeta_k[i]  *
      mean(dat_shells$B) ^ post_r$eta_b[i] 
  )
  psi <-   (mean(dat_shells$L)) ^ post_r$xi[i]
  R <- exp (  log(post_r$alpha[i] * phi * psi) + (post_r$sigma[i]^2 /2))
  samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_r$alpha[i] * phi * psi), post_r$sigma[i])
  points(jitter(seq(0,3,0.1)) * mean(dc_shellppl$age), samp_data, col = col.alpha("orange", 0.1), pch = 16)
  lines( seq(0,3,0.1) * mean(dc_shellppl$age),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}



curve(x^post_r$zeta_k[1], xlim = c(0,3), ylim = c(0, 2))
for(i in 1:350) curve(x^post_r$zeta_k[i], xlim = c(0,3), col = col.alpha("lightblue", 0.5),add = T)
curve(x^mean(post_r$zeta_k), xlim = c(0,3), col = "darkblue", add = T)
points(dat_shells$K, rep(1, length(dat_shells$K)), col = col.alpha("darkblue", 0.4), pch = 16)
curve(x^post_r$eta_b[1], xlim = c(0,3), ylim = c(0, 2))
for(i in 1:350) curve(x^post_r$eta_b[i], xlim = c(0,3), col = col.alpha("lightblue", 0.5),add = T)
curve(x^mean(post_r$eta_b), xlim = c(0,3), col = "darkblue", add = T)
points(dat_shells$B, rep(1, length(dat_shells$K)), col = col.alpha("darkblue", 0.4), pch = 16)


hist(max(dat_shells$K) ^ post_r$zeta_k - min(dat_shells$K) ^ post_r$zeta_k)
hist(max(dat_shells$B) ^ post_r$eta_b - min(dat_shells$B) ^ post_r$eta_b)
hist((max(dat_shells$K) ^ post_r$zeta_k - min(dat_shells$K) ^ post_r$zeta_k) - (max(dat_shells$B) ^ post_r$eta_b - min(dat_shells$B) ^ post_r$eta_b) )
hist((max(dat_shells$B) ^ post_r$eta_b - min(dat_shells$B) ^ post_r$eta_b) - (max(dat_shells$K) ^ post_r$zeta_k - min(dat_shells$K) ^ post_r$zeta_k))

sum((max(dat_shells$K) ^ post_r$zeta_k - min(dat_shells$K) ^ post_r$zeta_k) 
  - (max(dat_shells$B) ^ post_r$eta_b - min(dat_shells$B) ^ post_r$eta_b) 
  > 0)/length(post_r$eta_b)


par(mfrow = c(1,2))
plot(precis(m_ra ))
plot(precis(m_sa ))




dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$age),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]
dat_shells <- list(
     N = nrow(dc_shellppl),
     M = nrow(dc_shells),
     A = dc_shellppl$age[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$age),
     R = as.numeric(dc_shells$returns)/1000,
     L = dc_shells$lenght_min/mean(dc_shells$lenght_min),
     ID_i= as.integer(as.factor(as.character(dc_shells$anonymeID)))
   )
m_ro <- cstan( file= "models/Returns_o.stan" , data=dat_shells , chains=3, cores = 3 )
post_r <- extract.samples(m_ro)
plot(precis(m_ro))

plot(dat_shells$A, apply(post_r$phi, 2, mean), ylim=c(0,4))
for(i in 1:dat_shells$N) lines(rep(dat_shells$A[i],2), apply(post_r$phi, 2, PI)[,i])

curve((1-exp(- post_r$beta[1] * x  )) ^ post_r$gamma[1] , xlim = c(0,3), ylim = c(0,1), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_r$beta[i] * x  )) ^ post_r$gamma[i] , col = col.alpha("grey80", 0.5), add = T)

plot(dat_shells$A[dat_shells$ID_i], dat_shells$R, 
     xlim = c(0,3), ylim = c(0, max(dat$R)+1), pch = 16, col = col.alpha("grey40", 0.7),
     main = "gamma = lnorm(1,1)")
for(i in 1:150){
  phi <-  exp(apply(post_r$iota,1,mean )[i] ) * (
    (1-exp(- post_r$beta[i] * seq(0,3,0.1)  )) ^ post_r$gamma[i] 
  )
  psi <-   (mean(dat_shells$L)) ^ post_r$xi[i]
  R <-  exp(log(post_r$alpha[i] * phi * psi) + (post_r$sigma[i]^2 /2))
  samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_r$alpha[i] * phi * psi), post_r$sigma[i])
  points(jitter(seq(0,3,0.1)) * mean(dat_shells$A), samp_data, col = col.alpha("orange", 0.1), pch = 16)
  lines( seq(0,3,0.1) * mean(dat_shells$A),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}

###################################################################
###AAAAND IRT FOR KNOWLEDGE########################################
###################################################################
#simulated data
d<-sim_data(50,300, 200, zero = F, b_a = 1, g_a = 1, g_k = 1, g_b = 1)
dat <- list(
  N = d$N,
  M = d$M,
  Q = d$Q,
  A = d$A/mean(d$A),
  K = d$K/mean(d$K),
  Y = d$Y,
  B = d$B/mean(d$B),
  R = d$R,
  L = d$L/mean(d$L),
  ID_i = d$ID_ind
)
m_simirt <- cstan( file= "models/irt_only.stan" , data=dat , chains=3, cores = 3 )

#real data
m_rairt <- cstan( file= "models/Returns_allknow.stan" , data=dat_shells , chains=3, cores = 3 )
m_sairt <- cstan( file= "models/Success_allknow.stan" , data=dat_traps , chains=3, cores = 3 )
m_rhirt <- cstan( file= "models/Returns_irt_height.stan" , data=dat_shells , chains=3, cores = 3 )
m_rgirt <- cstan( file= "models/Returns_irt_grip.stan" , data=dat_shells , chains=3, cores = 3 )
m_shirt <- cstan( file= "models/Success_irt_height.stan" , data=dat_traps , chains=3, cores = 3 )
m_sgirt <- cstan( file= "models/Success_irt_grip.stan" , data=dat_traps , chains=3, cores = 3 )
par(mfrow = c(3,2), mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
models <- c("m_rairt", "m_sairt", "m_rhirt", "m_shirt", "m_rgirt", "m_sgirt")
for( i in 1:6) plot( precis( get (models [i])))

#shells
#actually to force stuff with priors seems to just work
post_irt <- extract.samples(m_rairt)
curve(inv_logit(mean(apply(post_irt$a, 2, mean)) * ( x - mean(apply(post_irt$b, 2, mean)))), 
      xlim = c(-6, 6), ylim = c(0, 1), 
      xlab = "latent scale", ylab = "p answer", 
      cex.lab=1.5, cex.axis=1.5,
      col = col.alpha("#1482ac", 1))
for (i in 1:dat_shells$Q) {
  curve(inv_logit(apply(post_irt$a, 2, mean)[i] * ( x - apply(post_irt$b, 2, mean)[i])), 
        col = col.alpha("#1482ac", 0.5),
        add = TRUE)}
points(apply(post_irt$K, 2, mean), jitter(rep(0.5, dat_shells$N)), 
       col = col.alpha("darkorange", 0.7), 
       pch = 16, cex=1.5)

#messy age effect
plot(dc_shellppl$age, apply(post_irt$K, 2, mean),
     xlab = "age", ylab = "irt-estimated knowledge", 
     col = col.alpha("darkorange", 0.7), pch = 16)
for (i in 1:ncol(post_irt$K)) {
  lines(rep(dc_shellppl$age[i], 2), apply(post_irt$K, 2, PI)[,i], col = col.alpha("grey80", 0.7) )  
}

#irt results are VERY VERY similar to simply counting number of items named
plot(dat_shells$K, apply(post_irt$K, 2, mean),
     xlim = c(0,2), ylim = c(-0.2,2),
     xlab = "knowledge based on number of items named", ylab = "irt-estimated knowledge", 
     col = col.alpha("darkorange", 0.7), pch = 16)
abline(0,1,col = col.alpha("#1482ac", 0.5) )
for (i in 1:ncol(post_irt$K)) {
     lines(rep(dat_shells$K[i], 2), apply(post_irt$K, 2, PI)[,i], col = col.alpha("grey80", 0.7) )  
}

par(mfrow = c(2,2))
plot(precis(m_rairt))
plot(dat_shells$A, apply(post_irt$phi, 2, mean), ylim=c(0,4))
for(i in 1:dat_shells$N) lines(rep(dat_shells$A[i],2), apply(post_irt$phi, 2, PI)[,i])

curve((1-exp(- post_irt$beta[1] * x  )) ^ post_irt$gamma[1] , xlim = c(0,3), ylim = c(0,2), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_irt$beta[i] * x  )) ^ post_irt$gamma[i] , col = col.alpha("grey80", 0.5), add = T)
for (i in 1:150) curve((1-exp(- post_irt$beta[i] * x  )) ^ post_irt$gamma[i]  *
                         mean(dat_shells$K) ^ post_irt$zeta_k[i], col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_irt$beta[i] * x  )) ^ post_irt$gamma[i]  *
                         mean(dat_shells$K) ^ post_irt$zeta_k[i] , col = col.alpha(shellcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_irt$beta[i] * x  )) ^ post_irt$gamma[i]  *
                         min(dat_shells$K) ^ post_irt$zeta_k[i], col = col.alpha(trapcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_irt$beta[i] * x  )) ^ post_irt$gamma[i]  *
                         max(dat_shells$K) ^ post_irt$zeta_k[i] , col = col.alpha(trapcol, 0.1), add = T)

plot(dc_shellppl$age[dat_shells$ID_i], dat_shells$R, 
     xlim = c(0,25), ylim = c(0, max(dat_shells$R)+1), pch = 16, col = col.alpha("grey40", 0.7))
for(i in 1:150){
  phi <-  exp(apply(post_irt$iota,1,mean )[i] ) * (
    (1-exp(- post_irt$beta[i] * seq(0,3,0.1)  )) ^ post_irt$gamma[i] *
      mean(dat_shells$K) ^ post_irt$zeta_k[i]
  )
  psi <-   (mean(dat_shells$L)) ^ post_irt$xi[i]
  R <- exp (  log(post_irt$alpha[i] * phi * psi) + (post_irt$sigma[i]^2 /2))
  samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_irt$alpha[i] * phi * psi), post_irt$sigma[i])
  points(jitter(seq(0,3,0.1)) * mean(dc_shellppl$age), samp_data, col = col.alpha("orange", 0.1), pch = 16)
  lines( seq(0,3,0.1) * mean(dc_shellppl$age),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}
par(mfrow = c(1,1))


#traps
post_s <- extract.samples(m_sairt)
curve(inv_logit(mean(apply(post_s$a, 2, mean)) * ( x - mean(apply(post_s$b, 2, mean)))), 
      xlim = c(-6, 6), ylim = c(0, 1), 
      xlab = "latent scale", ylab = "p answer", 
      cex.lab=1.5, cex.axis=1.5,
      col = col.alpha("#1482ac", 1))
for (i in 1:dat_traps$Q) {
  curve(inv_logit(apply(post_s$a, 2, mean)[i] * ( x - apply(post_s$b, 2, mean)[i])), 
        col = col.alpha("#1482ac", 0.5),
        add = TRUE)}
points(apply(post_s$K, 2, mean), jitter(rep(0.5, dat_traps$N)), 
       col = col.alpha("darkorange", 0.7), 
       pch = 16, cex=1.5)


#irt results are VERY VERY similar to simply counting number of items named
plot(dat_traps$K, apply(post_s$K, 2, mean),
     xlim = c(0,2), ylim = c(-0.2,2),
     xlab = "knowledge based on number of items named", ylab = "irt-estimated knowledge", 
     col = col.alpha("darkorange", 0.7), pch = 16)
abline(0,1,col = col.alpha("#1482ac", 0.5) )
for (i in 1:ncol(post_s$K)) {
     lines(rep(dat_traps$K[i], 2), apply(post_s$K, 2, PI)[,i], col = col.alpha("grey80", 0.7) )  
}

par(mfrow = c(2,2))
plot(precis(m_sairt))
plot(dat_traps$A, apply(post_s$phi, 2, mean), ylim=c(0,4))
for(i in 1:dat_traps$N) lines(rep(dat_traps$A[i],2), apply(post_s$phi, 2, PI)[,i])

curve((1-exp(- post_s$beta[1] * x  )) ^ post_s$gamma[1] , xlim = c(0,3), ylim = c(0,2), col = col.alpha("grey80", 0.5))
for (i in 1:150) curve((1-exp(- post_s$beta[i] * x  )) ^ post_s$gamma[i] , col = col.alpha("grey80", 0.5), add = T)
for (i in 1:150) curve((1-exp(- post_s$beta[i] * x  )) ^ post_s$gamma[i]  *
                         min(post_s$K[i,]) ^ post_s$zeta_k[i], col = col.alpha(trapcol, 0.1), add = T)
for (i in 1:150) curve((1-exp(- post_s$beta[i] * x  )) ^ post_s$gamma[i]  *
                         max(post_s$K[i,]) ^ post_s$zeta_k[i] , col = col.alpha(trapcol, 0.1), add = T)

plot(dc_trapppl$age[dat_traps$ID_i], dat_traps$S, 
     xlim = c(0,25), ylim = c(0, 1), pch = 16, col = col.alpha("grey40", 0.7))
for(i in 1:150){
  phi <-  exp(apply(post_s$iota,1,mean )[i] ) * (
    (1-exp(- post_s$beta[i] * seq(0,3,0.1)  )) ^ post_s$gamma[i] *
      mean(apply(post_s$K, 1, mean)) ^ post_s$zeta_k[i]
  )
  psi <-   (mean(dat_traps$L)) ^ post_s$xi[i]
  p <- 1 - exp ( - post_s$alpha[i] * phi * psi)
  samp_data <- rbern(length(seq(0,3,0.1)),  p)
  points(jitter(seq(0,3,0.1)) * mean(dc_trapppl$age), samp_data, col = col.alpha("orange", 0.1), pch = 16)
  lines( seq(0,3,0.1) * mean(dc_trapppl$age),  p, col = col.alpha(trapcol, 0.1), lwd = 2)
}
par(mfrow = c(1,1))

#########################
#tides###################
#########################
#donow5 check effects of the two with counterfactuals and duration
#with meters
m_tide_m <- cstan( file= "models/Returns_irt_tide.stan" , data=dat_shells , chains=3, cores = 3 )

#with integral
dat_shells$tide <- dc_shells$tide_integral
m_tide_i <- cstan( file= "models/Returns_irt_tide.stan" , data=dat_shells , chains=3, cores = 3 )


precis(m_tide_m)
precis(m_tide_i)

########################################################
#trying to make these cobb duglas behave as cobb douglas
########################################################

#drafts
m_shell_1k <- cstan( file= "models/2_shell_1k.stan" , data=dat_shells , chains=3, cores = 3 )
m_trap_1k <- cstan( file= "models/2_trap_1k.stan" , data=dat_traps , chains=3, cores = 3 )

dat_shells$H <- 1 + dat_shells$H
dat_traps$H <- 1 + dat_traps$H
m_shell_1h <- cstan( file= "models/2_shell_all.stan" , data=dat_shells , chains=3, cores = 3 )
m_trap_1h <- cstan( file= "models/2_trap_all.stan" , data=dat_traps , chains=3, cores = 3 )

dat_shells$H <- dc_shellppl$height / mean(dc_shellppl$height)
dat_traps$H <- dc_trapppl$height / mean(dc_trapppl$height)
dat_shells$G <- 1 + dat_shells$G
dat_traps$G <- 1 + dat_traps$G
m_shell_1g <- cstan( file= "models/2_shell_all.stan" , data=dat_shells , chains=3, cores = 3 )
m_trap_1g <- cstan( file= "models/2_trap_all.stan" , data=dat_traps , chains=3, cores = 3 )

par(mfrow = c(4,2), mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
models <- c("m_shell_all", "m_trap_all", "m_shell_1k", "m_trap_1k", "m_shell_1h", "m_trap_1h", "m_shell_1g", "m_trap_1g")
for( i in 1:8) plot( precis( get (models [i])))

for( j in 1:8){
  post <- extract.samples(get (models [j]))
  if(str_detect(models [j], "shell")){
      plot(dc_shellppl$age[dat_shells$ID_i], dat_shells$R, 
       xlim = c(0,age_plot), ylim = c(0, max(dat_shells$R)+1), pch = 16, col = col.alpha("grey40", 0.7))
      for(i in 1:150){
        phi <-  exp(apply(post$iota,1,mean )[i] ) * (
          (1-exp(- post$beta[i] * seq(0,3,0.1)  )) ^ post$gamma[i] *
            mean(post$K[i,]) ^ post$zeta_k[i] *
            1 ^ post$eta_h[i] *
            1 ^ post$theta_g[i]
        )
        psi <-   (mean(dat_shells$L)) ^ post$xi[i] * exp(post$tau[i] * mean(dat_shells$tide))
        R <- exp (  log(post$alpha[i] * phi * psi) + (post$sigma[i]^2 /2))
        samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post$alpha[i] * phi * psi), post$sigma[i])
        points(jitter(seq(0,3,0.1)) * mean(dc_shellppl$age), samp_data, col = col.alpha("orange", 0.1), pch = 16)
        lines( seq(0,3,0.1) * mean(dc_shellppl$age),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
    }
  }else{
      plot(dc_trapppl$age[dat_traps$ID_i], dat_traps$S, 
        xlim = c(0,age_plot), ylim = c(0, 1), pch = 16, col = col.alpha("grey40", 0.7))
      for(i in 1:150){
        phi <-  exp(apply(post$iota,1,mean )[i] ) * (
          (1-exp(- post$beta[i] * seq(0,3,0.1)  )) ^ post$gamma[i] *
            mean(post$K[i,]) ^ post$zeta_k[i]*
            1 ^ post$eta_h[i] *
            1 ^ post$theta_g[i]
        )
        psi <-   (mean(dat_traps$L)) ^ post$xi[i]
        p <- 1 - exp ( - post$alpha[i] * phi * psi)
        samp_data <- rbern(length(seq(0,3,0.1)),  p)
        points(jitter(seq(0,3,0.1)) * mean(dc_trapppl$age), samp_data, col = col.alpha("orange", 0.1), pch = 16)
        lines( seq(0,3,0.1) * mean(dc_trapppl$age),  p, col = col.alpha(trapcol, 0.1), lwd = 2)
      }

    }
}



# Maybe garbage?
# #fit to data
# plot(dat_shells$A[dat_shells$ID_ind] * mean(dc_shellppl$age), dat_shells$R, 
#      xlab = "Age", ylab = "Kg shells", main = titles[j],
#      xlim = c(0,30), pch = 16, col = col.alpha("grey40", 0.7))
# for(i in 1:150){
# phi <-  exp(apply(post_r$id_v,1,mean )[i]) * (
#    (1-exp(- post_r$beta[i] * seq(0,3,0.1)  )) ^ post_r$gamma[i] *
#    ifelse (length(post_r$zeta_k) > 0, (mean(dat_shells$K) ) ^ post_r$zeta_k[i] , 1 ) *
#    ifelse (length(post_r$eta_b) > 0, (mean(dat_shells$B) ) ^ post_r$eta_b[i] , 1 )
#   )
# psi <-   (mean(dat_shells$L)) ^ post_r$lambda[i]
# R <- exp (  log(post_r$alpha[i] * phi * psi) + (post_r$sigma[i]^2 /2))
# samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_r$alpha[i] * phi * psi), post_r$sigma[i])
# points(jitter(seq(0,3,0.1)) * mean(dc_shellppl$age), samp_data, col = col.alpha("orange", 0.1), pch = 16)
# lines( seq(0,3,0.1) * mean(dc_shellppl$age),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
# }
# 
# plot(jitter(dat_traps$A[dat_traps$ID_ind]) * mean(dc_trapppl$age), jitter(dat_traps$S, factor = 0.1), 
#      xlab = "Age", ylab = "Prob trap success", main = titles[j],
#      xlim = c(0,40), pch = 16, col = col.alpha("grey40", 0.2))
# for(i in 1:150){
# phi <-  exp(apply(post_s$id_v,1,mean )[i]) * (
#    (1-exp(- post_s$beta[i] * seq(0,3,0.1) ) ) ^ post_s$gamma[i] *
#    ifelse (length(post_s$zeta_k) > 0, (mean(dat_traps$K) ) ^ post_s$zeta_k[i] , 1 ) *
#    ifelse (length(post_s$eta_b) > 0, (mean(dat_traps$B) ) ^ post_s$eta_b[i] , 1 )
#   )
# psi <-   mean(dat_traps$L) ^ post_s$lambda[i]
# p <- 1 - exp  ( - post_s$alpha[i] * phi * psi )
# samp_data <- rbern(length(seq(0,3,0.1)), 1 - exp  ( -post_s$alpha[i]  * phi * psi))
# points(jitter(seq(0,3,0.1)) * mean(dc_trapppl$age), samp_data, col = col.alpha("lightgreen", 0.1), pch = 16)
# lines( seq(0,3,0.1) * mean(dc_trapppl$age),  p, col = col.alpha(trapcol, 0.1), lwd = 2)
# }

# par(mfrow = c(1,3))
# 
# 
# m_rak <- cstan( file= "models/Returns_allknow.stan" , data=dat_shells , chains=3, cores = 3 )




#######################################################################################
#REPRODUCE EHBEA PLOTS
#######################################################################################
quant <- c(0.9, 0.6, 0.3)
#age effect
m_ro <- cstan( file= "models/Returns_o.stan" , data=dat_shells , chains=3, cores = 3 )
m_so <- cstan( file= "models/Success_o.stan" , data=dat_traps , chains=3, cores = 3 )
m_rk <- cstan( file= "models/Returns_k.stan" , data=dat_shells , chains=3, cores = 3 )
m_sk <- cstan( file= "models/Success_k.stan" , data=dat_traps , chains=3, cores = 3 )
m_rb <- cstan( file= "models/Returns_b.stan" , data=dat_shells , chains=3, cores = 3 )
m_sb <- cstan( file= "models/Success_b.stan" , data=dat_traps , chains=3, cores = 3 )

post_r <- extract.samples(m_ro)
post_s <- extract.samples(m_so)
png("../plots/age_only.png", height = 6, width = 10, units = "in", res = 500)
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,30), ylim = c(0,1), 
         xlab = "Age", ylab = "Proportion improvement")#
    phi <- (1-exp(-median(post_r$beta) * seq_trait )) ^ median(post_r$gamma)
    lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
    phi <- (1-exp(-median(post_s$beta) * seq_trait )) ^ median(post_s$gamma)
    lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)
    for (q in 1:3){
          mu_phi <- sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta * x )) ^ post_r$gamma, quant[q]))
          shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.1))
          mu_phi <- sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta * x )) ^ post_s$gamma, quant[q]))
          shade(mu_phi, seq_trait * mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))
    }#q
dev.off()

#png("../plots/age_only.png", height = 6, width = 10, units = "in", res = 500)
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,30), ylim = c(0,1), 
     xlab = "Age", ylab = "Proportion improvement")#
phi <- (1-exp(-median(post_r$beta) * seq_trait )) ^ median(post_r$gamma)
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta) * seq_trait )) ^ median(post_s$gamma)
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)
for (q in 1:3){
  mu_phi <- sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta * x )) ^ post_r$gamma, quant[q]))
  shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.1))
  mu_phi <- sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta * x )) ^ post_s$gamma, quant[q]))
  shade(mu_phi, seq_trait * mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))
}#q
#dev.off()

#by min and max data for knowledge 
post_r <- extract.samples(m_rk)
post_s <- extract.samples(m_sk)
png("../plots/age_knowledge_2.png", height = 6, width = 10, units = "in", res = 500)
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,30), ylim = c(0,1.3), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta) * seq_trait  )) ^ median(post_r$gamma) *
       (mean(dat_shells$K)  ) ^ median(post_r$zeta_k) 
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta) * seq_trait  )) ^ median(post_s$gamma) *
       ( mean(dat_traps$K)  ) ^ median(post_s$zeta_k) 
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_r$beta) * seq_trait  )) ^ median(post_r$gamma) *
       ( min(dat_shells$K) ) ^ median(post_r$zeta_k)
mu_phi[2,] <- (1-exp(-median(post_r$beta) * seq_trait  )) ^ median(post_r$gamma) *
       (max(dat_shells$K)  ) ^ median(post_r$zeta_k)
shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.2))

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_s$beta) * seq_trait  )) ^ median(post_s$gamma) *
       ( min(dat_shells$K)  ) ^ median(post_s$zeta_k)
mu_phi[2,] <- (1-exp(-median(post_s$beta) * seq_trait  )) ^ median(post_s$gamma) *
       ( max(dat_shells$K)  ) ^ median(post_s$zeta_k)
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.2))
dev.off()

#by min and max data for body
post_r <- extract.samples(m_rb)
post_s <- extract.samples(m_sb)
png("../plots/age_height_2.png", height = 6, width = 10, units = "in", res = 500)
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,30), ylim = c(0,1.3), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta) * seq_trait  )) ^ median(post_r$gamma) *
       (mean(dat_shells$B)  ) ^ median(post_r$eta_b) 
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta) * seq_trait  )) ^ median(post_s$gamma) *
       ( mean(dat_traps$B)  ) ^ median(post_s$eta_b) 
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_r$beta) * seq_trait  )) ^ median(post_r$gamma) *
       ( min(dat_shells$B) ) ^ median(post_r$eta_b)
mu_phi[2,] <- (1-exp(-median(post_r$beta) * seq_trait  )) ^ median(post_r$gamma) *
       (max(dat_shells$B)  ) ^ median(post_r$eta_b)
shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.2))

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_s$beta) * seq_trait  )) ^ median(post_s$gamma) *
       ( min(dat_shells$B)  ) ^ median(post_s$eta_b)
mu_phi[2,] <- (1-exp(-median(post_s$beta) * seq_trait  )) ^ median(post_s$gamma) *
       ( max(dat_shells$B)  ) ^ median(post_s$eta_b)
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.2))
dev.off()
