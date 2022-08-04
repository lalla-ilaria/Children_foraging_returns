library(rethinking)
library(dagitty)
library(rlist)
library(tidyverse)
real_data <- list.load("2_data_preparation/processed_data.RData")
tide_data <- read.csv("2_data_preparation/tide_data.csv")

#####################
#PRIOR PLOTS
#####################
nsamp <-100
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

png("../plots/prior_age_only.png", height = 3, width = 4, units = "in", res = 500)
#plot prior predictive simulation
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,3), ylim = c (0,1), 
     xlab = "Age", ylab = "Proportion improvement",
     xaxt="n")
axis(1, at=seq(0,3,0.5),labels=seq(0,30,5))
#calculate per sample
for(i in 1:nsamp){
  phi[i,] <- (1-exp(-beta[i] * AGE  )) ^ gamma[i]
  lines( AGE,  phi[i,], col = col.alpha("cornflowerblue", 0.3))
}
lines (AGE, (1-exp(-1 * AGE  )) ^ 1,
       col = col.alpha("darkblue", 0.7))
text(1.5, 0.8,expression(paste(beta,"=1,",gamma,"=1",sep = "")))
lines (AGE, (1-exp(-1 * AGE  )) ^ 10,
       col = col.alpha("darkblue", 0.7))
text(2.5, 0.4,expression(paste(beta,"=1,",gamma,"=10",sep = "")))
lines (AGE, (1-exp(-3 * AGE  )) ^ 1,
       col = col.alpha("darkblue", 0.7))
text(0.4, 0.7,expression(paste(beta,"=3,",gamma,"=1",sep = "")))
dev.off()


png("../plots/prior_age_and_trait_4.png", height = 3, width = 4, units = "in", res = 500)
#plot prior predictive simulation
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,3), ylim = c (0,2), 
     xlab = "Age", ylab = "Proportion improvement",
     xaxt="n")
axis(1, at=seq(0,3,0.5),labels=seq(0,30,5))
#calculate per sample
lines (AGE, (1-exp(-1 * AGE  )) ^ 1,
       col = col.alpha("darkblue", 0.7))
text(0.2, 1.9,expression(paste(beta,"=1,",gamma,"=1",sep = "")))
lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( 0.7 ) ^ 1,
       col = col.alpha("cornflowerblue", 0.7))
text(2.7, 0.65,expression(paste(zeta,"=1",sep = "")))

lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( 1.3 ) ^ 1,
       col = col.alpha("cornflowerblue", 0.7))

lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( 0.7 ) ^ 2,
       col = col.alpha("lightblue", 0.7))
lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( 1.3 ) ^ 2,
       col = col.alpha("lightblue", 0.7))
text(2.7, 0.45,expression(paste(zeta,"=2",sep = "")))

lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( 0.7 ) ^ 0.3,
       col = col.alpha("dodgerblue4", 0.7))
lines (AGE, (1-exp(-1 * AGE  )) ^ 1 * 
         ( 1.3 ) ^ 0.3,
       col = col.alpha("dodgerblue4", 0.7))
text(2.7, 0.83,expression(paste(zeta,"=0.3",sep = "")))
dev.off()


for(i in 1:nsamp){
  phi[i,] <- (1-exp(-beta[i] * AGE  )) ^ gamma[i] *
    ( mean(K) ) ^ zeta_k[i]
  lines( AGE,  phi[i,], col = col.alpha("lightblue", 0.3))
}




###############################################
#QUADRYPTYCS###################################
###############################################
#Add these kinds of plots, check data, clean and make pretty

##########################################################################
#PREPARE REAL DATA
##########################################################################
#Real data
dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$knowledge),]
dc_shellppl <- dc_shellppl[complete.cases(dc_shellppl$height),]
dc_shell_k <- real_data$shell_k[which(rownames(real_data$shell_k) %in% dc_shellppl$anonymeID),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]
dat_shells <- list(
  N = nrow(dc_shellppl),
  M = nrow(dc_shells),
  Q = ncol(dc_shell_k),
  A = dc_shellppl$age[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$age),
  Y = dc_shell_k,
  K = dc_shellppl$knowledge[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$knowledge),
  B = dc_shellppl$height[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$height),
  H = dc_shellppl$height[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$height),
  G = dc_shellppl$grip[order(dc_shellppl$anonymeID)] / mean(dc_shellppl$grip),
  R = as.numeric(dc_shells$returns)/1000,
  L = dc_shells$lenght_min/mean(dc_shells$lenght_min),
  ID_i= as.integer(as.factor(as.character(dc_shells$anonymeID)))
)

dc_trapppl <- real_data$trap_ppl[complete.cases(real_data$trap_ppl$knowledge),]
dc_traps <- real_data$traps[which(real_data$traps$anonymeID %in% dc_trapppl$anonymeID),]
dc_trap_k <- real_data$trap_k[which(rownames(real_data$trap_k) %in% dc_trapppl$anonymeID),]
dat_traps <- list(
  N = nrow(dc_trapppl),
  M = nrow(dc_traps),
  Q = ncol(dc_trap_k),
  A = dc_trapppl$age[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$age),
  Y = dc_trap_k,
  K = dc_trapppl$knowledge[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$knowledge),
  B = dc_trapppl$height[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$height),
  H = dc_trapppl$height[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$height),
  G = dc_trapppl$grip[order(dc_trapppl$anonymeID)] / mean(dc_trapppl$grip),
  S = as.numeric(dc_traps$success),
  L = dc_traps$lenght_hour/mean(dc_traps$lenght_hour),
  ID_i= as.integer(as.factor(as.character(dc_traps$anonymeID)))
)



#mix cobb douglas

m_ra <- cstan( file= "models/Returns_all.stan" , data=dat_shells , chains=3, cores = 3 )
# m_rab <- cstan( file= "models/Returns_allbody.stan" , data=dat_shells , chains=3, cores = 3 )
m_sa <- cstan( file= "models/Success_all.stan" , data=dat_traps , chains=3, cores = 3 )

par(mfrow = c(2,2))
post_s <- extract.samples(m_sa)
post_r <- extract.samples(m_ra)
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

##########################
#TIDE PLOTS
##########################
png("../plots/tide_distance_kernel.png", height = 10, width = 12, units = "cm", res = 500)
  par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
  curve(exp(-0.07 * (x) ^2 ) *-1  , 
        xlim = c(-6,6), ylim = c(-1,0),
        xlab = "hours from max low tide", ylab = "proportion of high tide",
        col = "#1482ac", lwd = 2)
  points(c(-3.1, 3.1), c(-0.5,-0.5), pch = 16, cex = 1.5, col = "orange")
  text(0, -0.5, label = "mid-tide")
  arrows(-1, -0.5, -2.9, -0.5, length = 0.1)
  arrows(1, -0.5, 2.9, -0.5, length = 0.1)
dev.off()


png("../plots/tide_height_trip.png", height = 10, width = 12, units = "cm", res = 500)
  par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
  high_tide <- 3 #in meters. Consider searching for high tide value for each day #todo
  
  #plot tide depth and range of tide heights when foraging
  curve((exp(-0.07 * (x) ^2 ) * - 4 ) + high_tide , 
        xlim = c(-6,6), ylim = c(-1,4),
        xlab = "hours from max low tide", ylab = "depth of tide (m)",
        col = "white")
  #add shade for each time range of foraging trips
  for(i in 1:nrow(tide_data)){
    polygon(as.numeric( c( rep (tide_data$tide_start[i], 2), rep( tide_data$tide_end[i], 2))), 
            c(-1,4,4,-1), 
            col = col.alpha("lightblue", 0.1), border = NA)
  }
  #draw line of zero height for reference (low astronomical mean water level)
  abline(h = 0, col = "grey80" )
  #draw tide for each foraging trip
  for(i in 1:nrow(tide_data)){
    curve((exp(-0.07 * (x) ^2 ) * (tide_data$tide_height[i] - high_tide)) + high_tide, col = col.alpha("#1482ac", 0.3), add = TRUE)
  }
dev.off()

png("../plots/tide_height&avg.png", height = 10, width = 12, units = "cm", res = 500)
  ggplot(tide_data, aes(x = tide_height, y = avg_tide_depth, col = trip_length)) +
  geom_point( )+
  labs(x = "max tide height", y = "average tide height", col = "Duration (min)")+
  theme_classic()
dev.off()