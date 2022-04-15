library(rethinking)
library(rlist)
source("1_simulation/1b_simulation.R")
real_data <- list.load("2_data_preparation/processed_data.RData")
trapcol <- "#5ca81e"
shellcol <-  "#ea5914"
seq_trait <- seq(0,3,0.01)

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
beta_a <- rlnorm(nsamp, 0, 1)#prior
gamma_a <- rlnorm(nsamp, 0, 1)#prior
#sigma lognormal
sigma <- rexp(nsamp, 1)
#ADD KNOWLEDGE
beta_ak <- rlnorm(nsamp, 0, 1)#prior
gamma_k <- rlnorm(nsamp, 0, 1)#prior
K <- matrix(NA, nsamp, length(AGE)) 
for(i in 1:nsamp){
  K[i,] <- beta_ak[i]*AGE #trait
}
beta_ab <- rlnorm(nsamp, 0, 1)#prior
gamma_b <- rlnorm(nsamp, 0, 1)#prior
B <- matrix(NA, nsamp, length(AGE)) 
for(i in 1:nsamp){
  B[i,] <- beta_ab[i]*AGE #trait
}

#plot prior predictive simulation
plot(NULL, xlim = c(0,3), ylim = c (0,3), 
     xlab = "Age", ylab = "Proportion improvement")
#calculate per sample
for(i in 1:nsamp){
  phi[i,] <- (1-exp(-beta_a[i] * AGE  )) ^ gamma_a[i] 
  lines( AGE,  phi[i,], col = col.alpha("cornflowerblue", 0.3))
}
for(i in 1:nsamp){
  phi[i,] <- (1-exp(-beta_a[i] * AGE  )) ^ gamma_a[i] *
    ( mean(K) ) ^ gamma_k[i]
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
d<-sim_data(100,300, zero = F, b_a = 1, b_k = 0.1, b_b = 0.1)
d<-sim_data(100,300, zero = F, b_k = 0.5, b_b = 0.1)
plot(d$A, d$K)
plot(d$K[d$ID_ind], d$R)
dat <- list(
  N = d$N,
  M = d$M,
  A = d$A/mean(d$A),
  K = d$K/mean(d$K),
  B = d$B/mean(d$B),
  R = d$R/500,
  L = d$L/mean(d$L),
  ID_ind = d$ID_ind
)
m_ras <- cstan( file= "models/Returns_all4.stan" , data=dat , chains=3, cores = 3 )
post_r <- extract.samples(m_ras)

plot(d$A[d$ID_ind], d$R, 
     xlim = c(0,20), ylim = c(0, max(dat$R)), pch = 16, col = col.alpha("grey40", 0.7))
for(i in 1:150){
phi <-  apply(post_r$id_v,1,mean )[i] + (
   log(1-exp(- (post_r$beta_a[i] + 
                  post_r$beta_k[i] * mean(dat$K) +
                  post_r$beta_b[i] * mean(dat$B))
             * seq(0,3,0.1)  )) * post_r$gamma_a[i] 
  )
psi <-   log(mean(dat$L)) * post_r$lambda[i]
R <- exp (  post_r$alpha[i] + phi + psi + (post_r$sigma[i]^2 /2))
samp_data <- rlnorm(length(seq(0,3,0.1)),  post_r$alpha[i]  + phi + psi, post_r$sigma[i])
points(jitter(seq(0,3,0.1)) * mean(d$A), samp_data, col = col.alpha("orange", 0.1), pch = 16)
lines( seq(0,3,0.1) * mean(d$A),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}

plot(NULL, xlim = c(0,3), ylim = c(0,1))
for(i in 1:150){
  phi <- (1-exp(-post_r$beta_a[i] * seq_trait  )) ^ post_r$gamma_a[i]
  lines( seq_trait ,  phi, col = col.alpha("cornflowerblue", 0.4))
}
for(i in 1:150){
  phi <- (1-exp(-( post_r$beta_a[i] + post_r$beta_k[i] * mean(dat$K))  * seq_trait  )) ^ post_r$gamma_a[i]
  lines( seq_trait ,  phi, col = col.alpha("turquoise", 0.5))
}
for(i in 1:150){
  phi <- (1-exp(-( post_r$beta_a[i] + post_r$beta_b[i] * mean(dat$B))  * seq_trait  )) ^ post_r$gamma_a[i]
  lines( seq_trait ,  phi, col = col.alpha("lightblue", 0.5))
}


d<-sim_data(100,300, zero = T, b_a = 0.2, g_ak = 4)
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
  phi <- (1-exp(-post_s$beta_a[i] * seq_trait  )) ^ post_s$gamma_a[i]
  lines( seq_trait ,  phi, col = col.alpha("cornflowerblue", 0.4))
}
for(i in 1:150){
  phi <- (1-exp(-post_s$beta_k[i] * seq_trait  )) ^ post_s$gamma_k[i]
  lines( seq_trait ,  phi, col = col.alpha("turquoise", 0.5))
}
for(i in 1:150){
  phi <- (1-exp(-post_s$beta_b[i] * seq_trait  )) ^ post_s$gamma_b[i]
  lines( seq_trait ,  phi, col = col.alpha("lightblue", 0.5))
}




#mix cobb douglas

m_ra <- cstan( file= "models/Returns_all.stan" , data=dat_shells , chains=3, cores = 3 )
m_sa <- cstan( file= "models/Success_all.stan" , data=dat_traps , chains=3, cores = 3 )


par(mfrow = c(1,2))
plot(precis(m_ra ))
plot(precis(m_sa ))


post_r <- extract.samples(m_ra)
post_s <- extract.samples(m_sa)

#fit to data
plot(dat_shells$A[dat_shells$ID_ind] * mean(dc_shellppl$age), dat_shells$R, 
     xlab = "Age", ylab = "Kg shells", main = titles[j],
     xlim = c(0,30), pch = 16, col = col.alpha("grey40", 0.7))
for(i in 1:150){
phi <-  exp(apply(post_r$id_v,1,mean )[i]) * (
   (1-exp(- post_r$beta_a[i] * seq(0,3,0.1)  )) ^ post_r$gamma_a[i] *
   ifelse (length(post_r$beta_k) > 0, (mean(dat_shells$K) ) ^ post_r$gamma_k[i] , 1 ) *
   ifelse (length(post_r$beta_b) > 0, (mean(dat_shells$B) ) ^ post_r$gamma_b[i] , 1 )
  )
psi <-   (mean(dat_shells$L)) ^ post_r$lambda[i]
R <- exp (  log(post_r$alpha[i] * phi * psi) + (post_r$sigma[i]^2 /2))
samp_data <- rlnorm(length(seq(0,3,0.1)),  log(post_r$alpha[i] * phi * psi), post_r$sigma[i])
points(jitter(seq(0,3,0.1)) * mean(dc_shellppl$age), samp_data, col = col.alpha("orange", 0.1), pch = 16)
lines( seq(0,3,0.1) * mean(dc_shellppl$age),  R, col = col.alpha(shellcol, 0.1), lwd = 2)
}

plot(jitter(dat_traps$A[dat_traps$ID_ind]) * mean(dc_trapppl$age), jitter(dat_traps$S, factor = 0.1), 
     xlab = "Age", ylab = "Prob trap success", main = titles[j],
     xlim = c(0,40), pch = 16, col = col.alpha("grey40", 0.2))
for(i in 1:150){
phi <-  exp(apply(post_s$id_v,1,mean )[i]) * (
   (1-exp(- post_s$beta_a[i] * seq(0,3,0.1) ) ) ^ post_s$gamma_a[i] *
   ifelse (length(post_s$beta_k) > 0, (mean(dat_traps$K) ) ^ post_s$gamma_k[i] , 1 ) *
   ifelse (length(post_s$beta_b) > 0, (mean(dat_traps$B) ) ^ post_s$gamma_b[i] , 1 )
  )
psi <-   mean(dat_traps$L) ^ post_s$lambda[i]
p <- 1 - exp  ( - post_s$alpha[i] * phi * psi )
samp_data <- rbern(length(seq(0,3,0.1)), 1 - exp  ( -post_s$alpha[i]  * phi * psi))
points(jitter(seq(0,3,0.1)) * mean(dc_trapppl$age), samp_data, col = col.alpha("lightgreen", 0.1), pch = 16)
lines( seq(0,3,0.1) * mean(dc_trapppl$age),  p, col = col.alpha(trapcol, 0.1), lwd = 2)
}

par(mfrow = c(1,3))

#######################################################################################
#REPRODUCE EHBEA PLOTS
#######################################################################################
#age effect
plot(NULL, xlim = c(0,30), ylim = c(0,1), 
         xlab = "Age", ylab = "Proportion improvement", main = titles[j])#
    phi <- (1-exp(-median(post_r$beta_a) * seq_trait )) ^ median(post_r$gamma_a)
    lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
    phi <- (1-exp(-median(post_s$beta_a) * seq_trait )) ^ median(post_s$gamma_a)
    lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)
    for (q in 1:3){
          mu_phi <- sapply ( seq_trait , function (x) PI ((1-exp(-post_r$beta_a * x )) ^ post_r$gamma_a, quant[q]))
          shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(shellcol, 0.1))
          mu_phi <- sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta_a * x )) ^ post_s$gamma_a, quant[q]))
          shade(mu_phi, seq_trait * mean(dc_trapppl$age), col = col.alpha(trapcol, 0.1))
    }#q

#by min and max data for knowledge 
plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (mean(dat_shells$K)  ) ^ median(post_r$gamma_k) 
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       ( mean(dat_traps$K)  ) ^ median(post_s$gamma_k) 
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       ( min(dat_shells$K) ) ^ median(post_r$gamma_k)
mu_phi[2,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (max(dat_shells$K)  ) ^ median(post_r$gamma_k)
shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.2))

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       ( min(dat_shells$K)  ) ^ median(post_s$gamma_k)
mu_phi[2,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       ( max(dat_shells$K)  ) ^ median(post_s$gamma_k)
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.2))

#by min and max data for body
plot(NULL, xlim = c(0,30), ylim = c(0,1), xlab = "Age", ylab = "Proportion improvement")
phi <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (mean(dat_shells$B)  ) ^ median(post_r$gamma_b) 
lines( seq_trait * mean(dc_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
phi <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       ( mean(dat_traps$B)  ) ^ median(post_s$gamma_b) 
lines( seq_trait * mean(dc_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       ( min(dat_shells$B) ) ^ median(post_r$gamma_b)
mu_phi[2,] <- (1-exp(-median(post_r$beta_a) * seq_trait  )) ^ median(post_r$gamma_a) *
       (max(dat_shells$B)  ) ^ median(post_r$gamma_b)
shade(mu_phi, seq_trait* mean(dc_shellppl$age), col = col.alpha(shellcol, 0.2))

mu_phi <- matrix(nrow = 2, ncol = length(seq_trait)) 
mu_phi[1,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       ( min(dat_shells$B)  ) ^ median(post_s$gamma_b)
mu_phi[2,] <- (1-exp(-median(post_s$beta_a) * seq_trait  )) ^ median(post_s$gamma_a) *
       ( max(dat_shells$B)  ) ^ median(post_s$gamma_b)
shade(mu_phi, seq_trait* mean(dc_trapppl$age), col = col.alpha(trapcol, 0.2))

