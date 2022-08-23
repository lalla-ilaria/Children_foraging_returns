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

####################################################################################################################
#FROM 1_ANALYSIS DRAFTS OF DATA IMPUTATION
####################################################################################################################
##########################################################################
#PREPARE DATA -  input knowledge - draft
##########################################################################
#SHELLS
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),                       #n individuals in total sample
  M = nrow(d_shells),                         #n trip/person
  ID_i= d_shells$index_id,                    #index of person of trip 
  has_foraging = ifelse(d_shellppl$data == "shells", 1, 0),
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = (d_shellppl$age / mean(d_shellppl$age)),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #knowledge data
  has_knowledge = ifelse(is.na(d_shellppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  Q = ncol(d_shell_k),                        #n items in freelist
  answers = d_shell_k                       #all answers from freelist
)
dat_shells[["answers"]][is.na(dat_shells[["answers"]])] <- -999

m_shell_knowledge<- cstan( file= "models/3_shells_imputk.stan" , data=dat_shells , 
                           chains=3, cores = 3, iter = 500 )
#TRAPS
dat_traps <- list(
  #foraging data
  N = nrow(d_trapppl),                       #n individuals in total sample
  M = nrow(d_traps),                         #n trip/person
  ID_i= d_traps$index_id,                    #index of person of trip 
  has_foraging = ifelse(d_trapppl$data == "traps", 1, 0),
  success = d_traps$success,                 #whether trap captured something
  age = (d_trapppl$age / mean(d_trapppl$age)),
  sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_traps$lenght_hour/mean(d_traps$lenght_hour),
  #knowledge data
  has_knowledge = ifelse(is.na(d_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  Q = ncol(d_trap_k),                        #n items in freelist
  answers = d_trap_k                       #all answers from freelist
)

dat_traps[["answers"]][is.na(dat_traps[["answers"]])] <- -999


m_traps_knowedge <- cstan( file= "models/3_traps_imputk.stan" , data=dat_traps , 
                           chains=3, cores = 3, iter = 500 )


##########################################################################
#PREPARE DATA -  input height - draft
##########################################################################
#NB need not to estimate knowledge for the people we have only anthropometrics
#but make one list of data 
#SHELLS
#create data frame
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells

#remove people for whom we have no height nor foraging returns
d_shellppl <- d_shellppl[-which(is.na(d_shellppl$height)&
                                  d_shellppl$data != "shells"),]
#removing people for whom we have no age
d_shellppl <- d_shellppl[-which(is.na(d_shellppl$age)&
                                  d_shellppl$data != "shells"),]

#add index variables
#index and sort all individuals so we can loop across them
d_shellppl$index_id <- as.integer(as.factor(d_shellppl$anonymeID))
d_shellppl <- d_shellppl[order(d_shellppl$index_id),]
#add index for individuals in the foraging data
for ( i in 1:nrow(d_shells)){
  d_shells$index_id[i] <- d_shellppl$index_id[which ( d_shellppl$anonymeID == d_shells$anonymeID[i])]
}

dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),                       #n individuals in total sample
  M = nrow(d_shells),                         #n trip/person
  ID_i= d_shells$index_id,                    #index of person of trip 
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = (d_shellppl$age / mean(d_shellppl$age)),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #height data
  has_height = ifelse(is.na(d_shellppl$height), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  height = d_shellppl$height/mean(d_shellppl$height, na.rm = TRUE),
  min_height = 50/mean(d_shellppl$height, na.rm = TRUE)
)
dat_shells[["height"]][is.na(dat_shells[["height"]])] <- -999

m_shells_height <- cstan( file= "models/4_shell_height.stan" , data=dat_shells , 
                          chains=1, cores = 1, iter = 500 )

#TRAPS
#create data frame
d_trapppl <- real_data$trap_ppl
d_traps <- real_data$traps

#remove people for whom we have no height nor foraging returns
d_trapppl <- d_trapppl[-which(is.na(d_trapppl$height)&
                                d_trapppl$data != "traps"),]
#removing people for whom we have no age
d_trapppl <- d_trapppl[-which(is.na(d_trapppl$age)&
                                d_trapppl$data != "traps"),]

#add index variables
#index and sort all individuals so we can loop across them
d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
#add index for individuals in the foraging data
for ( i in 1:nrow(d_traps)){
  d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
}

dat_traps <- list(
  #foraging data
  N = nrow(d_trapppl),                       #n individuals in total sample
  M = nrow(d_traps),                         #n trip/person
  ID_i= d_traps$index_id,                    #index of person of trip 
  success = d_traps$success,                 #whether trap captured something
  age = (d_trapppl$age / mean(d_trapppl$age)),
  sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_traps$lenght_hour/mean(d_traps$lenght_hour),
  #height data
  has_height = ifelse(is.na(d_trapppl$height), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  height = d_trapppl$height/mean(d_trapppl$height, na.rm = TRUE),
  min_height = 50/mean(d_trapppl$height, na.rm = TRUE)
)
dat_traps[["height"]][is.na(dat_traps[["height"]])] <- -999

m_traps_height <- cstan( file= "models/4_trap_height.stan" , data=dat_traps , 
                         chains=1, cores = 1, iter = 500 )


##########################################################################
#PREPARE DATA -  input height AND knowledge - draft
##########################################################################
#NB need not to estimate knowledge for the people we have only anthropometrics
#but make one list of data 
#SHELLS
#create data frame
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

dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),                       #n individuals in total sample
  M = nrow(d_shells),                         #n trip/person
  ID_i= d_shells$index_id,                    #index of person of trip 
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = (d_shellppl$age / mean(d_shellppl$age)),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #height data
  has_height = ifelse(is.na(d_shellppl$height), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  height = d_shellppl$height/mean(d_shellppl$height, na.rm = TRUE),
  min_height = 50/mean(d_shellppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
  #knowledge data
  has_knowledge = ifelse(is.na(d_shellppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  Q = ncol(d_shell_k),                        #n items in freelist
  answers = d_shell_k                         #all answers from freelist
)
dat_shells[["answers"]][is.na(dat_shells[["answers"]])] <- -999
dat_shells[["height"]][is.na(dat_shells[["height"]])] <- -999

m_shells_height <- cstan( file= "models/5_shell_heightandknowledge.stan" , data=dat_shells , 
                          chains=1, cores = 1, iter = 500 )


##########################################################################
#PREPARE DATA -  input grip - draft
##########################################################################
#NB need not to estimate knowledge for the people we have only anthropometrics
#but make one list of data 
#SHELLS
#create data frame
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells

#remove people for whom we have no height nor foraging returns
d_shellppl <- d_shellppl[-which(is.na(d_shellppl$grip)&
                                  d_shellppl$data != "shells"),]

#add index variables
#index and sort all individuals so we can loop across them
d_shellppl$index_id <- as.integer(as.factor(d_shellppl$anonymeID))
d_shellppl <- d_shellppl[order(d_shellppl$index_id),]
#add index for individuals in the foraging data
for ( i in 1:nrow(d_shells)){
  d_shells$index_id[i] <- d_shellppl$index_id[which ( d_shellppl$anonymeID == d_shells$anonymeID[i])]
}

dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),                       #n individuals in total sample
  M = nrow(d_shells),                         #n trip/person
  ID_i= d_shells$index_id,                    #index of person of trip 
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = (d_shellppl$age / mean(d_shellppl$age)),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #grip data
  has_grip = ifelse(is.na(d_shellppl$grip), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  grip = d_shellppl$grip/mean(d_shellppl$grip, na.rm = TRUE),
)
dat_shells[["grip"]][is.na(dat_shells[["grip"]])] <- -999

m_shells_grip <- cstan( file= "models/6_shell_grip.stan" , data=dat_shells , 
                        chains=1, cores = 1, iter = 500 )

######################################
#other stuff -scroll through and check
######################################

precis(m_draft)
post <- extract.samples(m_traps1)

#make plots
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_knowledge == "1", "deepskyblue4", "darkorange3")

#check standardized  knowledge
plot(NULL, xlim = c(1,dat_shells$N), ylim = c(-10, 4), 
     xlab = "index", ylab = "standardized knowledge")
for (i in 1:dat_shells$N) {
  points ( rep(i, 750 ),
           post$knowledge[,i] ,
           col = rep(ifelse(dat_shells$has_knowledge[i] == 1, col.alpha("deepskyblue", 0.3), col.alpha("darkgoldenrod1", 0.3))))
}
points(1:dat_shells$N, apply(post$knowledge, 2, mean), 
       ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_shells$N) {
  lines ( rep(i, 2 ),
          PI(post$knowledge[,i]) ,
          col = rep(ifelse(dat_shells$has_knowledge[i] == 1, "deepskyblue4", "darkorange3")))
}

#check knowledge by age
plot(d_shellppl$knowledge, apply(post$knowledge, 2, mean))
year_eff <- apply(post$delta_j[,], 1, cumsum)
year_seq <- seq(0,4,0.2)
plot(x = dat_shells$age, 
     y = apply(post$knowledge, 2, mean), 
     xlim = c(0,4),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,1] * ( 1 - exp(-post$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,2] * ( 1 - exp(-post$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

#plot of returns
phi <-  mean(post$iota) +
  mean(post$gamma) * log(1-exp(- mean(post$beta) * seq_trait  )) +
  mean(post$zeta_k)* mean(post$knowledge) 
psi <-  mean(post$xi) * mean(log(dat_shells$duration)) +
  mean(post$tau)* mean(dat_shells$tide) 
samp_data <- rlnorm(length(seq_trait),  
                    mean(log(post$alpha)) + phi + psi, 
                    mean(post$sigma))

plot(jitter(seq_trait) * mean(d_shellppl$age), samp_data, 
     xlim = c(0,age_plot), ylim = c(0, max(dat_shells$returns)+1), 
     xlab = "Age", ylab = "kg shellfish",
     pch = 16, col = col.alpha("orange", 0.2))
for(i in 1:150){
  phi <-  apply(post$iota,1,mean )[i] +
    post$gamma[i] * log(1-exp(- post$beta[i] * seq_trait  )) +
    post$zeta_k[i]*apply(post$knowledge, 1, mean)[i] 
  
  psi <-  post$xi[i] * mean(log(dat_shells$duration)) + 
    post$tau[i] * mean(dat_shells$tide) 
  R <- exp (  log(post$alpha[i]) + phi + psi + 
                (post$sigma[i]^2 /2))
  lines( seq_trait * mean(d_shellppl$age),  R, 
         col = col.alpha(shellcol, 0.2), lwd = 1)
}
points(jitter(d_shellppl$age[d_shells$index_id], amount = 0.4), 
       as.numeric(d_shells$returns)/1000,
       pch = 16, cex = 0.8, col = col.alpha(othercol, 0.4))

#diff knowledge outcome scale
post <- extract.samples(m_draft4)
phi <-  apply(post$iota,1,mean )  +
  post$gamma * log(1-exp(- post$beta * 20  )) + 
  post$zeta_k* min(apply(post$knowledge, 2, mean))   

psi <-  post$xi * mean(log(dat_shells$duration)) + 
  post$tau* mean(dat_shells$tide)
min_k <- exp (mean(log(post$alpha)) + phi + psi + 
                (mean(post$sigma)^2 /2))
phi <-  apply(post$iota,1,mean )  +
  post$gamma * log(1-exp(- post$beta * 20  )) + 
  post$zeta_k* max(apply(post$knowledge, 2, mean))   

psi <-  post$xi * mean(log(dat_shells$duration)) + 
  post$tau* mean(dat_shells$tide)
max_k <- exp (mean(log(post$alpha)) + phi + psi + 
                (mean(post$sigma)^2 /2))
diff_k1 <- max_k - min_k
plot(density(diff_k1))



#TRAPS
d_trapppl <- real_data$trap_ppl
d_traps <- real_data$traps
d_trap_k <- real_data$trap_k

#remove people for whom we have only anthropometric data
d_trapppl <- d_trapppl %>% filter(!data == "anthropometrics")
d_trapppl$age[which(is.na(d_trapppl$age))] <- 45
d_trap_k <- d_trap_k[1:nrow(d_trapppl),]

#add index variables
#index and sort all individuals so we can loop across them
d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
for ( i in 1:nrow(d_traps)){
  d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
}
#sort knowledge data
d_trap_k <- d_trap_k[ order(row.names(d_trap_k)), ]

dat_traps <- list(
  #foraging data
  N = nrow(d_trapppl),                       #n individuals in total sample
  M = nrow(d_traps),                         #n trip/person
  ID_i= d_traps$index_id,                    #index of person of trip 
  has_foraging = ifelse(d_trapppl$data == "traps", 1, 0),
  has_knowledge = ifelse(is.na(d_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  success = d_traps$success,                 #whether trap captured something
  age = (d_trapppl$age / mean(d_trapppl$age)),
  sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_traps$lenght_hour/mean(d_traps$lenght_hour),
  #knowledge data
  Q = ncol(d_trap_k),                        #n items in freelist
  answers = d_trap_k                       #all answers from freelist
)

dat_traps[["answers"]][is.na(dat_traps[["answers"]])] <- -999


m_traps_knowedge <- cstan( file= "models/3_trap.stan" , data=dat_traps , 
                           chains=3, cores = 3, iter = 500 )

precis(m_traps_knowedge)
post <- extract.samples(m_traps_knowedge)
sex_col <- ifelse(dat_traps$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_traps$has_knowledge == "1", "deepskyblue4", "darkorange3")

#check standardized  knowledge
plot(NULL, xlim = c(1,dat_traps$N), ylim = c(-10, 4), 
     xlab = "index", ylab = "standardized knowledge")
for (i in 1:dat_traps$N) {
  points ( rep(i, 750 ),
           post$knowledge[,i] ,
           col = rep(ifelse(dat_traps$has_knowledge[i] == 1, col.alpha("deepskyblue", 0.3), col.alpha("darkgoldenrod1", 0.3))))
}
points(1:dat_traps$N, apply(post$knowledge, 2, mean), 
       ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_traps$N) {
  lines ( rep(i, 2 ),
          PI(post$knowledge[,i]) ,
          col = rep(ifelse(dat_traps$has_knowledge[i] == 1, "deepskyblue4", "darkorange3")))
}

#check knowledge by age
plot(d_trapppl$knowledge, apply(post$knowledge, 2, mean))
year_seq <- seq(0,4,0.2)
plot(x = dat_traps$age, 
     y = apply(post$knowledge, 2, mean), 
     xlim = c(0,4),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_traps$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,1] * ( 1 - exp(-post$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,2] * ( 1 - exp(-post$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}


#plot of returns
phi <-  mean(post$iota) +
  mean(post$gamma) * log(1-exp(- mean(post$beta) * seq_trait  )) +
  mean(post$zeta_k)* mean(post$knowledge) 
psi <-  mean(post$xi) * mean(log(dat_shells$duration)) 
samp_data <- rbern(length(seq_trait),  
                   1 - exp (-mean(post$alpha) * exp(phi) * exp(psi)))

plot(jitter(seq_trait) * mean(d_trapppl$age), samp_data, 
     xlim = c(0,age_plot), ylim = c(0, 1), 
     xlab = "Age", ylab = "p capture",
     pch = 16, col = col.alpha("lawngreen", 0.2))
for(i in 1:150){
  phi <-  apply(post$iota,1,mean )[i] +
    post$gamma[i] * log(1-exp(- post$beta[i] * seq_trait  )) +
    post$zeta_k[i]*apply(post$knowledge, 1, mean)[i] 
  
  psi <-  post$xi[i] * mean(log(dat_traps$duration)) 
  p <- 1 - exp (-post$alpha[i] * exp(phi) * exp(psi) )
  lines( seq_trait * mean(d_trapppl$age),  p, 
         col = col.alpha(trapcol, 0.2), lwd = 1)
}
points(jitter(d_trapppl$age[d_traps$index_id], amount = 0.5), 
       jitter(d_traps$success, amount = 0.02), 
       pch = 16, cex = ifelse(d_traps$success == 1, 0.8, 0.6), 
       col = col.alpha(othercol, ifelse(d_traps$success == 1, 0.4, 0.1)))






################3
#################
#################

#SHELLS
#create data frame
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells
d_shell_k <- real_data$shell_k
d_trapppl <- real_data$trap_ppl
d_traps <- real_data$traps
d_trap_k <- real_data$trap_k


#remove people for whom we have only anthropometric data
d_shellppl <- d_shellppl %>% filter(!data == "anthropometrics")
#d_shellppl <- d_shellppl %>% filter(!data == "knowledge")
d_shell_k <- d_shell_k[1:nrow(d_shellppl),]
d_trapppl <- d_trapppl %>% filter(!data == "anthropometrics")
d_trapppl$age[which(is.na(d_trapppl$age))] <- 45
#d_trapppl <- d_trapppl %>% filter(!data == "knowledge")
d_trap_k <- d_trap_k[1:nrow(d_trapppl),]



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

#same for traps
d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
for ( i in 1:nrow(d_traps)){
  d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
}
d_trap_k <- d_trap_k[ order(row.names(d_trap_k)), ]

dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  has_foraging = ifelse(d_shellppl$data == "shells", 1, 0),#[1:40],
  has_knowledge = ifelse(is.na(d_shellppl$knowledge), 0, 1),#[1:40], #vector of 0/1 for whether knowledge has to be imputed
  returns = as.numeric(d_shells$returns)/1000,
  age = (d_shellppl$age / mean(d_shellppl$age)),#[1:40],
  age_int = d_shellppl$age,#[1:40],
  sex = ifelse(d_shellppl$sex == "m", 1, 2),#[1:40], #make vector of sexes 1 = male 2 = female
  #height 
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #knowledge data
  Q = 100,#ncol(d_shell_k),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  answers = d_shell_k[, 1:100], #all answers from freelist
  prior_dirichlet = rep( 0.5, length (1:56 ) -1 ) #prior for Dirichlet
)
for (i in 1:length(dat_shells)) dat_shells[[i]][is.na(dat_shells[[i]])] <- -999

#try tightening priors
#making a positive version of knowledge
#connect to foraging stuff
m_shell_irt1 <- cstan( file= "models/other_models/2_shell_k3.stan" , data=dat_shells , 
                       chains=3, cores = 3, iter = 500 )

#m <- cmdstan_model("models/2_shell_k2.stan" )
#standardizing knowledge IN the models???
post <- extract.samples(m_shell_irt)
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_knowledge == "1", "deepskyblue4", "darkorange3")
plot(NULL, xlim = c(1,dat_shells$N), ylim = c(0, 4))
for (i in 1:dat_shells$N) {
  points ( rep(i, 750 ),
           post$knowledge_st[,i] ,
           col = rep(ifelse(dat_shells$has_knowledge[i] == 1, col.alpha("deepskyblue", 0.3), col.alpha("darkgoldenrod1", 0.3))))
}
points(1:dat_shells$N, apply(post$knowledge_st, 2, mean), 
       ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_shells$N) {
  lines ( rep(i, 2 ),
          PI(post$knowledge_st[,i]) ,
          col = rep(ifelse(dat_shells$has_knowledge[i] == 1, "deepskyblue4", "darkorange3")))
}

plot(d_shellppl$knowledge, apply(post$knowledge_st, 2, mean))
year_eff <- apply(post$delta_j[,], 1, cumsum)
year_seq <- seq(0,4,0.2)
plot(x = dat_shells$age, 
     y = apply(post$knowledge, 2, mean), 
     xlim = c(0,4),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,1] * ( 1 - exp(-post$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,2] * ( 1 - exp(-post$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}
#plot of returns
post_s <- post
phi <-  exp(mean(post_s$iota) ) * (
  (1-exp(- mean(post_s$beta) * seq_trait  )) ^ mean(post_s$gamma))
psi <-   (mean(dat_shells$duration)) ^ mean(post_s$xi) * 
  exp(mean(dat_shells$tide) * mean(post_s$tau))
R <- exp (  log(mean(post_s$alpha) * phi * psi) + 
              (mean(post_s$sigma)^2 /2))
samp_data <- rlnorm(length(seq_trait),  
                    log(mean(post_s$alpha) * phi * psi), 
                    mean(post_s$sigma))

plot(jitter(seq_trait) * mean(d_shellppl$age), samp_data, 
     xlim = c(0,age_plot), ylim = c(0, max(dat_shells$returns)+1), 
     xlab = "Age", ylab = "kg shellfish",
     pch = 16, col = col.alpha("orange", 0.2))
for(i in 1:150){
  phi <-  exp(apply(post_s$iota,1,mean )[i] ) * (
    (1-exp(- post_s$beta[i] * seq_trait  )) ^ post_s$gamma[i]
  )
  psi <-   (mean(dat_shells$duration)) ^ post_s$xi[i] * 
    exp(mean(dat_shells$tide) * post_s$tau[i])
  R <- exp (  log(post_s$alpha[i] * phi * psi) + 
                (post_s$sigma[i]^2 /2))
  lines( seq_trait * mean(d_shellppl$age),  R, 
         col = col.alpha(shellcol, 0.2), lwd = 1)
}
points(jitter(d_shellppl$age[d_shells$index_id], amount = 0.4), 
       as.numeric(d_shells$returns)/1000,
       pch = 16, cex = 0.8, col = col.alpha(othercol, 0.4))

#m_shell_all <- cstan( file= "models/2_shell_knowledgeinput.stan" , data=dat_shells , chains=3, cores = 3 )

#TRAPS
dat_traps <- list(
  #foraging data
  N = nrow(d_trapppl),                       #n individuals in total sample
  M = nrow(d_traps),                         #n trip/person
  ID_i= d_traps$index_id,                    #index of person of trip 
  has_foraging = ifelse(d_trapppl$data == "traps", 1, 0),
  has_knowledge = ifelse(is.na(d_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  success = d_traps$success,                 #whether trap captured something
  age = (d_trapppl$age / mean(d_trapppl$age)),
  sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_traps$lenght_hour/mean(d_traps$lenght_hour),
  #knowledge data
  Q = ncol(d_trap_k),                        #n items in freelist
  answers = d_trap_k                       #all answers from freelist
)
dat_traps <- list(
  #foraging data
  N = nrow(d_trapppl),                       #n individuals in total sample
  M = nrow(d_traps),                         #n trip/person
  ID_i= d_traps$index_id,                    #index of person of trip 
  has_foraging = ifelse(d_trapppl$data == "traps", 1, 0),
  has_knowledge = ifelse(is.na(d_trapppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  success = d_traps$success,                 #whether trap captured something
  age = (d_trapppl$age / mean(d_trapppl$age)),
  sex = ifelse(d_trapppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_traps$lenght_hour/mean(d_traps$lenght_hour),
  #knowledge data
  Q = ncol(d_trap_k),                        #n items in freelist
  answers = d_trap_k                       #all answers from freelist
)

dat_traps[["answers"]][is.na(dat_traps[["answers"]])] <- -999


m_traps1 <- cstan( file= "models/3_trap.stan" , data=dat_traps , 
                   chains=3, cores = 3, iter = 500 )


dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  ID_k = d_shell_k_ppl$index_id, #make vector of length W where is reported index ID of individuals in foraging df
  returns = as.numeric(d_shells$returns)/1000,
  age = d_shellppl$age / mean(d_shellppl$age),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  #height
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  knowledge_impute = ifelse(is.na(d_shellppl$knowledge), 1, 0), #vector of 0/1 for whether knowledge has to be imputed
  #knowledge data
  W = nrow(d_shell_k_ppl), #n of individuals for whom we have knowledge
  Q = ncol(d_knowledge$Y_l),
  O = 57 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age, #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2),
  answers = d_knowledge$Y_l, #all answers from freelist
  prior_dirichlet = rep( 0.5, length (0:40 ) -1 ) #prior for Dirichlet
)
m_shell_all <- cstan( file= "models/2_shell_knowledgeinput.stan" , data=dat_shells , chains=3, cores = 3 )
# for (i in 1:nrow(d_shellppl)) {
#   if(d_shellppl$anonymeID[i] %in% rownames(d_shell_k)) {
#     d_shellppl$index_k[i] <- which (rownames(d_shell_k) == d_shellppl$anonymeID[i])
#   } else {
#     d_shellppl$index_k[i] <- NA
#   }
# }
# for (i in 1:nrow(d_trapppl)) {
#   if(d_trapppl$anonymeID[i] %in% rownames(d_trap_k)) {
#     d_trapppl$index_k[i] <- which (rownames(d_trap_k) == d_trapppl$anonymeID[i])
#   } else {
#     d_trapppl$index_k[i] <- NA
#   }
# }

dat_shells <- list(
  #knowledge data
  W = 20,#nrow(d_shell_k_ppl), #n of individuals for whom we have knowledge
  Q = 50,#ncol(d_knowledge$Y_l),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age[1:20], #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2)[1:20],
  answers = d_knowledge$Y_l[1:20, 1:50], #all answers from freelist
  prior_dirichlet = rep( 2, length (1:56 ) -1 ) #prior for Dirichlet
)
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  ID_k = ifelse( is.na(d_shellppl$position_k), -999, d_shellppl$position_k) , #make vector of length N where is reported position of individuals for whom we have knowledge in the knowledge list
  returns = as.numeric(d_shells$returns)/1000,
  age = d_shellppl$age / mean(d_shellppl$age),
  age_int = d_shellppl$age ,
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  #height 
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  knowledge_impute = ifelse(is.na(d_shellppl$knowledge), 1, 0), #vector of 0/1 for whether knowledge has to be imputed
  #knowledge data
  W = nrow(d_shell_k), #n of individuals for whom we have knowledge
  Q = ncol(d_shell_k),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age, #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2),
  answers = d_shell_k, #all answers from freelist
  prior_dirichlet = rep( 0.2, length (1:56 ) -1 ) #prior for Dirichlet
)
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  ID_k = ifelse( is.na(d_shellppl$position_k), -999, d_shellppl$position_k) , #make vector of length N where is reported position of individuals for whom we have knowledge in the knowledge list
  returns = as.numeric(d_shells$returns)/1000,
  age = d_shellppl$age / mean(d_shellppl$age),
  age_int = d_shellppl$age ,
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  #height 
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  knowledge_impute = ifelse(is.na(d_shellppl$knowledge), 1, 0), #vector of 0/1 for whether knowledge has to be imputed
  #knowledge data
  W = nrow(d_shell_k_ppl), #n of individuals for whom we have knowledge
  Q = ncol(d_knowledge$Y_l),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age, #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2),
  answers = d_knowledge$Y_l, #all answers from freelist
  prior_dirichlet = rep( 0.2, length (1:56 ) -1 ) #prior for Dirichlet
)

#m_shell_irt <- cstan( file= "models/irt_only_3.stan" , data=dat_shells , chains=3, cores = 3 )


for (i in 1:150) {
  lines(x = 1:nrow(year_eff),  
        y = post$omega[i] + post$ro_age[i,1] * year_eff[,i], 
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:150) {
  lines(x = 1:nrow(year_eff),  
        y = post$omega[i] + post$ro_age[i,2] * year_eff[,i], 
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}
lines( x = 1:nrow(year_eff),  
       y = mean(post$omega) + mean(post$ro_age[,1]) * apply(year_eff, 1, mean), 
       type = "l", 
       col = col.alpha( boycol, alpha = 1), lwd = 1.5 )
lines( x = 1:nrow(year_eff),  
       y = mean(post$omega) + mean(post$ro_age[,2]) * apply(year_eff, 1, mean), 
       type = "l", 
       col = col.alpha( girlcol, alpha = 1), lwd = 1.5 )

#for height
# #diff height outcome scale
# phi <-  apply(post$iota,1,mean )  +
#         post$gamma * log(1-exp(- post$beta * 20  )) + 
#         post$eta_h* min(dat_shells$height[which(dat_shells$height >= 0)])   
# 
# psi <-  post$xi * mean(log(dat_shells$duration)) + 
#         post$tau* mean(dat_shells$tide)
# min_h <- exp (mean(log(post$alpha)) + phi + psi + 
#               (mean(post$sigma)^2 /2))
# phi <-  apply(post$iota,1,mean )  +
#         post$gamma * log(1-exp(- post$beta * 20  )) + 
#         post$eta_h* max(dat_shells$height)   
# 
# psi <-  post$xi * mean(log(dat_shells$duration)) + 
#         post$tau* mean(dat_shells$tide)
# max_h <- exp (mean(log(post$alpha)) + phi + psi + 
#               (mean(post$sigma)^2 /2))
# diff_h1 <- max_h - min_h
# plot(density(diff_h1))
# 
# 

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
phi_min_k <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))

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

