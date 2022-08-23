library(rethinking)
library(dagitty)
library(rlist)
library(tidyverse)
real_data <- list.load("2_data_preparation/processed_data.RData")
tide_data <- read.csv("2_data_preparation/tide_data.csv")
boycol  <- rgb(114/255, 181/255, 40/255) #"navyblue"
girlcol <- rgb(208/255, 29/255, 157/255) #"red3"
trapcol <- "#52b502" #"#5ca81e"
shellcol <-  "#ea5914"
othercol <- "grey30"##1482ac"
seq_trait <- seq(0,3,0.001)
age_plot <- 40

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

############################################################
#CALCULATE PRIORS FOR EFFECTS OF KNOWLEDGE ETC
############################################################
#load fit to extract knowledge values and overlay posterior
post_s <- extract.samples(m_shells_neg)
post_t <- extract.samples(m_traps_neg)

n_samp <- nrow(diffs_phi)
iota <- rnorm(n_samp, 0, 1)
alpha <- log(exp(rnorm(n_samp, 0, 1))) + 0.0001
beta <- rexp(n_samp, 1)
gamma <- rexp(n_samp, 1)
zeta_k <- rnorm(n_samp, 0, 1)
eta_h <- rnorm(n_samp, 0, 1)
theta_g <- rnorm(n_samp, 0, 1)
xi <- rnorm(n_samp, 0, 1)
tau <- rnorm(n_samp, 0, 1)
sigma <- rexp(n_samp, 1)

diffs_out <- data.frame(matrix(nrow=n_samp, ncol = 9))
colnames(diffs_out) <- c("sk","sh","sg","sd","st","tk","th","tg","td")
diffs_phi <- data.frame(matrix(nrow=n_samp, ncol = 6))
colnames(diffs_phi) <- c("sk","tk","sh","th","sg","tg")
#####
#SHELLS_prior
#####
phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* min(apply(post_s$knowledge, 2, mean))
psi <-  xi * mean(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
min_k <- exp (log(alpha) + phi_min + psi +
                (sigma^2 /2))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* max(apply(post_s$knowledge, 2, mean))
psi <-  xi * mean(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
max_k <- exp (log(alpha) + phi_max + psi +
                (sigma^2 /2))
diffs_out$sk <- max_k - min_k
diffs_phi$sk <- phi_max - phi_min

phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* min(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  xi * mean(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
min_h <- exp (log(alpha) + phi_min + psi +
                (sigma^2 /2))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* max(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  xi * mean(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
max_h <- exp (log(alpha) + phi_max + psi +
                (sigma^2 /2))
diffs_out$sh <- max_h - min_h
diffs_phi$sh <- phi_max - phi_min

phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* min(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  xi * mean(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
min_g <- exp (log(alpha) + phi_min + psi +
                (sigma^2 /2))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* max(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
phi_max <-  xi * mean(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
max_g <- exp (log(alpha) + phi_max + psi +
                (sigma^2 /2))
diffs_out$sg <- max_g - min_g
diffs_phi$sg <- phi_max - phi_min

phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  xi * min(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
min_d <- exp (log(alpha) + phi_min + psi +
                (sigma^2 /2))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  xi * max(log(dat_shells$duration)) +
  tau* mean(dat_shells$tide)
max_d <- exp (log(alpha) + phi_max + psi +
                (sigma^2 /2))
diffs_out$sd <- max_d - min_d

phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  xi * mean(log(dat_shells$duration)) +
  tau* min(dat_shells$tide)
min_t <- exp (log(alpha) + phi_min + psi +
                (sigma^2 /2))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  xi * mean(log(dat_shells$duration)) +
  tau* max(dat_shells$tide)
max_t <- exp (log(alpha) + phi_max + psi +
                (sigma^2 /2))
diffs_out$st <- min_t - max_t
#####
diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out)),
                                  rep("tide", nrow(diffs_out))),
                    kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
)

diffs_shells_prior <- diffs 




#TRAPS_prior
######
phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* min(apply(post_t$knowledge, 2, mean))
psi <-  xi * mean(log(dat_traps$duration)) 
min_k <- 1 - exp ( - alpha * exp(phi_min) * exp(psi))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* max(apply(post_t$knowledge, 2, mean))
psi <-  xi * mean(log(dat_traps$duration)) 
max_k <- 1 - exp ( - alpha * exp(phi_max) * exp(psi))
diffs_out$tk <- (max_k - min_k)[1:length(alpha)]
diffs_phi$tk <- (phi_max - phi_min)[1:length(alpha)]


phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* min(log(dat_traps$height), na.rm = TRUE) +
  theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  xi * mean(log(dat_traps$duration)) 
min_h <- 1 - exp ( - alpha * exp(phi_min) * exp(psi))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* max(log(dat_traps$height), na.rm = TRUE) +
  theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  xi * mean(log(dat_traps$duration)) 
max_h <- 1 - exp ( - alpha * exp(phi_max) * exp(psi))
diffs_out$th <- (max_h - min_h)[1:length(alpha)]
diffs_phi$th <- (phi_max - phi_min)[1:length(alpha)]

phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  theta_g* min(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  xi * mean(log(dat_traps$duration)) 
min_g <- 1 - exp ( - alpha * exp(phi_min) * exp(psi))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  theta_g* max(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  xi * mean(log(dat_traps$duration)) 
max_g <- 1 - exp ( - alpha * exp(phi_max) * exp(psi))
diffs_out$tg <- (max_g - min_g)[1:length(alpha)]
diffs_phi$tg <- (phi_max - phi_min)[1:length(alpha)]

phi_min <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  xi * min(log(dat_traps$duration)) 
min_d <- 1 - exp ( - alpha * exp(phi_min) * exp(psi))
phi_max <-  iota  +
  gamma * log(1-exp(- beta * 20  )) +
  eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  xi * max(log(dat_traps$duration)) 
max_d <- 1 - exp ( - alpha * exp(phi_max) * exp(psi))
diffs_out$td <- (max_d - min_d )[1:length(alpha)]
#####

diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out))),
                    p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
)

diffs_traps_prior <- diffs 


v_diffs_phi_prior <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi)),
                                        rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi))),
                          p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
                          foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))


######
#CALCULATE DIFFERENCE BETWEEN MINIMUM AND MAXIMUM ON OUTCOME SCALE
diffs_out <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 9))
colnames(diffs_out) <- c("sk","sh","sg","sd","st","tk","th","tg","td")
diffs_phi <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 6))
colnames(diffs_phi) <- c("sk","tk","sh","th","sg","tg")
######
#SHELLS_post
######
phi_min <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* min(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
min_k <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
phi_max <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* max(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
max_k <- exp (log(post_s$alpha) + phi_max + psi +
                (post_s$sigma^2 /2))
diffs_out$sk <- max_k - min_k
diffs_phi$sk <- phi_max - phi_min

phi_min <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* min(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
min_h <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
phi_max <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* max(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
max_h <- exp (log(post_s$alpha) + phi_max + psi +
                (post_s$sigma^2 /2))
diffs_out$sh <- max_h - min_h
diffs_phi$sh <- phi_max - phi_min

phi_min <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* min(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
min_g <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
phi_max <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* max(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
phi_max <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
max_g <- exp (log(post_s$alpha) + phi_max + psi +
                (post_s$sigma^2 /2))
diffs_out$sg <- max_g - min_g
diffs_phi$sg <- phi_max - phi_min

phi_min <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * min(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
min_d <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
phi_max <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * max(log(dat_shells$duration)) +
  post_s$tau* mean(dat_shells$tide)
max_d <- exp (log(post_s$alpha) + phi_max + psi +
                (post_s$sigma^2 /2))
diffs_out$sd <- max_d - min_d

phi_min <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* min(dat_shells$tide)
min_t <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
phi_max <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(dat_shells$duration)) +
  post_s$tau* max(dat_shells$tide)
max_t <- exp (log(post_s$alpha) + phi_max + psi +
                (post_s$sigma^2 /2))
diffs_out$st <- min_t - max_t
#####
diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out)),
                                  rep("tide", nrow(diffs_out))),
                    kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
)

diffs_shells <- diffs 
#TRAPS-post
#####
phi_min <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* min(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * mean(log(dat_traps$duration)) 
min_k <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
phi_max <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* max(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * mean(log(dat_traps$duration)) 
max_k <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
diffs_out$tk <- (max_k - min_k)[1:length(post_s$alpha)]
diffs_phi$tk <- (phi_max - phi_min)[1:length(post_s$alpha)]


phi_min <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* min(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * mean(log(dat_traps$duration)) 
min_h <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
phi_max <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* max(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * mean(log(dat_traps$duration)) 
max_h <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
diffs_out$th <- (max_h - min_h)[1:length(post_s$alpha)]
diffs_phi$th <- (phi_max - phi_min)[1:length(post_s$alpha)]

phi_min <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* min(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * mean(log(dat_traps$duration)) 
min_g <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
phi_max <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* max(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * mean(log(dat_traps$duration)) 
max_g <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
diffs_out$tg <- (max_g - min_g)[1:length(post_s$alpha)]
diffs_phi$tg <- (phi_max - phi_min)[1:length(post_s$alpha)]

phi_min <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * min(log(dat_traps$duration)) 
min_d <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
phi_max <-  apply(post_t$iota,1,mean )  +
  post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
  post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
  post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
  post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
psi <-  post_t$xi * max(log(dat_traps$duration)) 
max_d <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
diffs_out$td <- (max_d - min_d )[1:length(post_s$alpha)]
#####

diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out))),
                    p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
)

diffs_traps <- diffs 


out_shells <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(data = diffs_shells_prior, aes(x = variable, y = kg_shellfish), fill = col.alpha("orange", 0.4),
              colour = "orange") +
  geom_violin(fill = col.alpha(shellcol, 0.4),
              colour = shellcol) +
  labs(y = "kg shellfish difference", x = "")+
  ylim(-4, 6)+
  theme_classic()

png("../plots/diffout_shell_violin_pior.png", height = 10, width = 10, units = "cm", res = 500)
out_shells
dev.off()

out_traps <- ggplot(diffs_traps, aes(x = variable, y = p_success)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(data = diffs_traps_prior, aes(x = variable, y = p_success), fill = col.alpha("lawngreen", 0.4),
              colour = "lawngreen") +
  geom_violin(fill = col.alpha(trapcol, 0.3),
              colour = trapcol) +
  labs(y = "p_success difference", x = "")+
  ylim(-1, 1)+
  theme_classic()

png("../plots/diffout_trap_violin_pior.png", height = 10, width = 10, units = "cm", res = 500)
out_traps
dev.off()



v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi)),
                                        rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi))),
                          p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
                          foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))




phi_trait_prior <- ggplot(v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type, color = foraging_type) )+ 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(data = v_diffs_phi_prior, aes(x = variable, y = p_success, fill = foraging_type, color = foraging_type), trim=FALSE)+
  scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c(  col.alpha("orange", 0.4), col.alpha("lawngreen", 0.4 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c("orange", "lawngreen"),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  new_scale_fill()+
  new_scale_color()+
  geom_violin(data = v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type, color = foraging_type), trim=FALSE)+
  scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c( col.alpha(shellcol, 0.4), col.alpha(trapcol, 0.4 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c(shellcol, trapcol),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  labs(x = "", y = "phi difference")+
  ylim(-10, 10)+
  theme_classic()+
  theme(legend.position="top")

png("../plots/phi_trait_violin_pior.png", height = 8, width = 14, units = "cm", res = 500)
phi_trait_prior
dev.off()


png("../plots/alldiffs_combined_prior.png", height = 14, width = 16, units = "cm", res = 500)
plots <- align_plots(phi_trait_prior, out_shells, align = 'v', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], out_traps, labels = c('B', 'C'), label_size = 12)

# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
dev.off()

##################################################
#AGE ONLY - traps with average 
##################################################
post_t <- extract.samples(m_trap_age)
png("../plots/age_only_trap_mean.png", height = 8, width = 8, units = "cm", res = 500)
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)

#traps 
phi <-  mean(post_t$iota) +
  mean(post_t$gamma) * log(1-exp(- mean(post_t$beta) * seq_trait  )) 
psi <-  mean(post_t$xi) * (mean(log (dat_traps$duration))) 
p <- 1 - exp ( - mean(post_t$alpha) * exp(phi) * exp(psi))
samp_data <- rbern(length(seq_trait),  p)

#NB making plot with 13212 only (one of best hunters) to make plot with optimal situation to raise curve off zero
actor <- d_trapppl$index_id[which(d_trapppl$anonymeID == 13212)]
plot(jitter(seq_trait) * mean(d_trapppl$age), samp_data, 
     xlab = "Age", ylab = "p capture",
     xlim = c(0,age_plot), ylim = c(0, 1), 
     pch = 16, col = col.alpha("lawngreen", 0.3))
# #with ideal conditions (best actor, shortest time)
# for(i in 1:150){
#   phi <-  post_t$iota[i,actor] + 
#     post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait)) 
#   psi <-  post_t$xi[i] * log(min(dat_traps$duration))  
#   p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
#   lines( seq_trait * mean(d_trapppl$age),  p, 
#          col = col.alpha(trapcol, 0.2), lwd = 1)
# }
#with average actor and average time 
for(i in 1:150){
  phi <-  apply(post_t$iota,1,mean )[i] +
          post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait))
  psi <-  post_t$xi[i] * mean(log(dat_traps$duration))
  p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
  lines( seq_trait * mean(d_trapppl$age),  p,
         col = col.alpha(trapcol, 0.2), lwd = 1)
}
points(jitter(d_trapppl$age[d_traps$index_id], amount = 0.5), 
       jitter(d_traps$success, amount = 0.02), 
       pch = 16, cex = ifelse(d_traps$success == 1, 0.8, 0.6), 
       col = col.alpha(othercol, ifelse(d_traps$success == 1, 0.4, 0.1)))
dev.off()

##########################################################################
#AGE ONLY -TIDE MODEL
##########################################################################

for(i in 1:nrow(d_shells)){
  d_shells$age[i] <- d_shellppl$age[which (d_shellppl$anonymeID == d_shells$anonymeID[i])]
}
for(i in 1:nrow(d_shells)){
  d_shells$sex[i] <- d_shellppl$sex[which (d_shellppl$anonymeID == d_shells$anonymeID[i])]
}
sex_col <- ifelse(d_shells$sex == "m", boycol, girlcol)

dat_shells <- list(
  M = nrow(d_shells),
  age = d_shells$age / mean(d_shells$age),
  tide = d_shells$tide_avg_depth
)

m_tide <- cstan(file = "models/tide_age.stan", data = dat_shells)
post_tide <- extract.samples (m_tide)

png("../plots/age_and_tide.png", height = 3, width = 4, units = "in", res = 500)
plot(d_shells$age, d_shells$tide_avg_depth, pch = 16, col = sex_col)
for ( i in 1:250) lines( seq(0,4, 0.1) * mean(d_shells$age), 
                         post_tide$alpha[i] + post_tide$beta[i] * seq(0,4, 0.1),
                         col = col.alpha("#1482ac", 0.3))
dev.off()
##################################################
#N ITEMS
##################################################
png("../plots/n_items_age.png", height = 3, width = 4, units = "in", res = 500)
plot(d_shells$age[which(d_shells$n_item_types > 0)], 
     d_shells$n_item_types[which(d_shells$n_item_types > 0)],
     xlim = c(5,age_plot ), ylim = c(0, 10),
     xlab = "Age", ylab = "n types of shellfish",
     pch = 16, col = "#1482ac"  )
dev.off()

##########################################################################
#KNOWLEDGE ESTIMATION PLOTS - draft
##########################################################################

#SHELLS
post_s <- extract.samples(m_shell_knowledge)
#make plots
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_knowledge == "1", "deepskyblue4", "darkorange3")

#check standardized  knowledge
plot(NULL, xlim = c(1,dat_shells$N), ylim = c(-10, 4), 
     xlab = "index", ylab = "irt knowledge")
for (i in 1:dat_shells$N) {
  points ( rep(i, nrow(post_s$knowledge) ),
    post_s$knowledge[,i] ,
    col = rep(ifelse(dat_shells$has_knowledge[i] == 1, col.alpha("deepskyblue", 0.3), col.alpha("darkgoldenrod1", 0.3))))
}
points(1:dat_shells$N, apply(post_s$knowledge, 2, mean), 
     ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_shells$N) {
  lines ( rep(i, 2 ),
    PI(post_s$knowledge[,i]) ,
    col = rep(ifelse(dat_shells$has_knowledge[i] == 1, "deepskyblue4", "darkorange3")))
}

#estimated knowledge vs n of items listed
plot(d_shellppl$knowledge, apply(post_s$knowledge, 2, mean), 
     xlab = "n items listed", ylab = "irt knowledge",
     pch = 19, col = presence_col)

#check knowledge by age
year_seq <- seq(0,4,0.2)
plot(x = dat_shells$age * mean(d_shellppl$age, na.rm = TRUE), 
     y = apply(post_s$knowledge, 2, median), 
     xlim = c(0,40),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:250) {
  lines(x = year_seq * mean(d_shellppl$age, na.rm = TRUE),  
        y = post_s$omega[i] + post_s$chi[i,1] * ( 1 - exp(-post_s$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:250) {
  lines(x = year_seq * mean(d_shellppl$age, na.rm = TRUE),  
        y = post_s$omega[i] + post_s$chi[i,2] * ( 1 - exp(-post_s$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

#plot of returns
  phi <-  mean(post_s$iota) +
          mean(post_s$gamma) * log(1-exp(- mean(post_s$beta) * seq_trait  )) +
          mean(post_s$zeta_k)* mean(post_s$knowledge) 
  psi <-  mean(post_s$xi) * mean(log(dat_shells$duration)) +
          mean(post_s$tau)* mean(dat_shells$tide) 
  samp_data <- rlnorm(length(seq_trait),  
                      mean(log(post_s$alpha)) + phi + psi, 
                      mean(post_s$sigma))
  
  plot(jitter(seq_trait) * mean(d_shellppl$age), samp_data, 
       xlim = c(0,age_plot), ylim = c(0, max(dat_shells$returns)+1), 
       xlab = "Age", ylab = "kg shellfish",
       pch = 16, col = col.alpha("orange", 0.2))
  for(i in 1:150){
    phi <-  apply(post_s$iota,1,mean )[i] +
            post_s$gamma[i] * log(1-exp(- post_s$beta[i] * seq_trait  )) +
            post_s$zeta_k[i]*apply(post_s$knowledge, 1, mean)[i] 
    
    psi <-  post_s$xi[i] * mean(log(dat_shells$duration)) + 
            post_s$tau[i] * mean(dat_shells$tide) 
    R <- exp (  log(post_s$alpha[i]) + phi + psi + 
                  (post_s$sigma[i]^2 /2))
    lines( seq_trait * mean(d_shellppl$age),  R, 
           col = col.alpha(shellcol, 0.2), lwd = 1)
  }
  points(jitter(d_shellppl$age[d_shells$index_id], amount = 0.4), 
         as.numeric(d_shells$returns)/1000,
         pch = 16, cex = 0.8, col = col.alpha(othercol, 0.4))
  
#plot with age variation only
  plot(NULL, xlim = c(0,age_plot), ylim = c(0,1), 
       xlab = "Age", ylab = "proportion max foraging")#, main = "Age only"
  phi <- exp(median(post_s$gamma) * log(1-exp(-median(post_s$beta) * seq_trait  ))) 
  lines( seq_trait * mean(d_shellppl$age),  phi, col = col.alpha(shellcol, 1), lwd = 2)
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_s$gamma * log(1-exp(-post_s$beta * x ))), 0.95) )
  shade(mu_phi, seq_trait* mean(d_trapppl$age), col = col.alpha(shellcol, 0.1))
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_s$gamma * log(1-exp(-post_s$beta * x ))), 0.5) )
  shade(mu_phi, seq_trait* mean(d_trapppl$age), col = col.alpha(shellcol, 0.15))
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta * x )) ^ post_s$gamma, 0.3) )
  shade(mu_phi, seq_trait* mean(d_trapppl$age), col = col.alpha(shellcol, 0.15))
  for(i in 1:30){
    phi <- exp(post_s$gamma[i] * log(1-exp(-post_s$beta[i] * seq_trait )))
    lines( seq_trait * mean(d_shellppl$age),  phi, col = col.alpha(shellcol, 0.3))
  }

  
#TRAPS
post_t <- extract.samples(m_traps_knowedge)
#make plots
sex_col <- ifelse(dat_traps$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_traps$has_knowledge == "1", "deepskyblue4", "darkorange3")

#check standardized  knowledge
plot(NULL, xlim = c(1,dat_traps$N), ylim = c(-10, 4), 
     xlab = "index", ylab = "standardized knowledge")
for (i in 1:dat_traps$N) {
  points ( rep(i, 750 ),
    post_t$knowledge[,i] ,
    col = rep(ifelse(dat_traps$has_knowledge[i] == 1, col.alpha("deepskyblue", 0.3), col.alpha("darkgoldenrod1", 0.3))))
}
points(1:dat_traps$N, apply(post_t$knowledge, 2, mean), 
     ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_traps$N) {
  lines ( rep(i, 2 ),
    PI(post_t$knowledge[,i]) ,
    col = rep(ifelse(dat_traps$has_knowledge[i] == 1, "deepskyblue4", "darkorange3")))
}

#estimated knowledge vs n items listed
plot(d_trapppl$knowledge, apply(post_t$knowledge, 2, mean), 
     xlab = "n items listed", ylab = "irt knowledge",
     pch = 19, col = presence_col)

#check knowledge by age
year_seq <- seq(0,4,0.2)
plot(x = dat_traps$age, 
     y = apply(post_t$knowledge, 2, mean), 
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
        y = post_t$omega[i] + post_t$chi[i,1] * ( 1 - exp(-post_t$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post_t$omega[i] + post_t$chi[i,2] * ( 1 - exp(-post_t$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}


#plot of returns
  phi <-  mean(post_t$iota) +
          mean(post_t$gamma) * log(1-exp(- mean(post_t$beta) * seq_trait  )) +
          mean(post_t$zeta_k)* mean(post_t$knowledge) 
  psi <-  mean(post_t$xi) * mean(log(dat_shells$duration)) 
  samp_data <- rbern(length(seq_trait),  
                      1 - exp (-mean(post_t$alpha) * exp(phi) * exp(psi)))
  
  actor <- d_trapppl$index_id[which(d_trapppl$anonymeID == 13212)]
  plot(jitter(seq_trait) * mean(d_trapppl$age), samp_data, 
       xlim = c(0,age_plot), ylim = c(0, 1), 
       xlab = "Age", ylab = "p capture",
       pch = 16, col = col.alpha("lawngreen", 0.2))
  #with ideal conditions (best actor, shortest time)
  for(i in 1:150){
    phi <-  post_t$iota[i,actor] + 
            post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait)) +
            post_t$zeta_k[i]*apply(post_t$knowledge, 1, mean)[i] 
    psi <-  post_t$xi[i] * log(min(dat_traps$duration))  
    p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
    lines( seq_trait * mean(d_trapppl$age),  p, 
           col = col.alpha(trapcol, 0.2), lwd = 1)
  }
  #with mean actor and time
  for(i in 1:150){
    phi <-  apply(post_t$iota,1,mean )[i] +
            post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait  )) +
            post_t$zeta_k[i]*apply(post_t$knowledge, 1, mean)[i] 
    
    psi <-  post_t$xi[i] * mean(log(dat_traps$duration)) 
    p <- 1 - exp (-post_t$alpha[i] * exp(phi) * exp(psi) )
    lines( seq_trait * mean(d_trapppl$age),  p, 
           col = col.alpha(trapcol, 0.2), lwd = 1)
  }
 points(jitter(d_trapppl$age[d_traps$index_id], amount = 0.5), 
         jitter(d_traps$success, amount = 0.02), 
    pch = 16, cex = ifelse(d_traps$success == 1, 0.8, 0.6), 
    col = col.alpha(othercol, ifelse(d_traps$success == 1, 0.4, 0.1)))

 #age variation only
   plot(NULL, xlim = c(0,age_plot), ylim = c(0,1), 
       xlab = "Age", ylab = "proportion max foraging")#, main = "Age only"
  phi <-  exp(median(post_t$gamma) * log(1-exp(-median(post_t$beta) * seq_trait  ))) 
  lines( seq_trait * mean(d_trapppl$age),  phi, col = col.alpha(trapcol, 1), lwd = 2)
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_t$gamma * log(1-exp(-post_t$beta * x ))), 0.95) )
  shade(mu_phi, seq_trait* mean(d_trapppl$age), col = col.alpha(trapcol, 0.1))
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_t$gamma * log(1-exp(-post_t$beta * x ))), 0.5) )
  shade(mu_phi, seq_trait* mean(d_trapppl$age), col = col.alpha(trapcol, 0.15))
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_t$gamma * log(1-exp(-post_t$beta * x ))), 0.3) )
  shade(mu_phi, seq_trait* mean(d_trapppl$age), col = col.alpha(trapcol, 0.15))
  for(i in 1:30){
    phi <-  exp(post_t$gamma[i] * log(1-exp(-post_t$beta[i] * seq_trait )))
    lines( seq_trait * mean(d_trapppl$age),  phi, col = col.alpha(trapcol, 0.3))
  }
  
##########################################################################
#HEIGHT ESTIMATION - draft
##########################################################################
#SHELLS
#create data frame
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells

#remove people for whom we have no height nor foraging returns
d_shellppl <- d_shellppl[-which(is.na(d_shellppl$height)&
                                  d_shellppl$data != "shells"),]
#removing people for whom we have no age
d_shellppl$age[which (d_shellppl$anonymeID == "12588")] <- 7
d_shellppl$age[which (d_shellppl$anonymeID == "252472")] <- 30
d_shellppl$age[which (d_shellppl$anonymeID == "84993")] <- 15

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
  has_height = ifelse(is.na(d_shellppl$height), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  height = d_shellppl$height/mean(d_shellppl$height, na.rm = TRUE),
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = (d_shellppl$age / mean(d_shellppl$age)),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth
)
dat_shells[["height"]][is.na(dat_shells[["height"]])] <- -999

post <- extract.samples(m_shells_height)

sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_height == "1", "deepskyblue4", "darkorange3")
year_seq <- seq(0,4,0.2)

plot(x = dat_shells$age, 
     y = apply(post_s$height_merged, 2, mean), 
     xlim = c(0,4),
     xlab = "Age" , 
     ylab = "Height",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_height == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:250) {
  lines(x = year_seq ,  
        y = dat_shells$min_height + post_s$kappa[i,1] * ( 1 - exp(-post_s$lambda[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:250) {
  lines(x = year_seq ,  
        y = dat_shells$min_height + post_s$kappa[i,2] * ( 1 - exp(-post_s$lambda[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

#plot of returns
  phi <-  mean(post$iota) +
          mean(post$gamma) * log(1-exp(- mean(post$beta) * seq_trait  )) +
          mean(post$eta_h)* mean(post$height_merged) 
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
            post$eta_h[i]*apply(post$height_merged, 1, mean)[i] 
    
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


  #TRAPS
  #create data frame
  
  post <- extract.samples(m_traps_height)
  
  sex_col <- ifelse(dat_traps$sex == "1", boycol, girlcol)
  presence_col <- ifelse(dat_traps$has_height == "1", "deepskyblue4", "darkorange3")
  year_seq <- seq(0,4,0.2)
  
  plot(x = dat_traps$age, 
       y = apply(post$height_merg, 2, mean), 
       xlim = c(0,4),
       xlab = "Age" , 
       ylab = "Height",
       cex.lab=1.8 , 
       cex.axis=1.8 ,
       pch = ifelse(dat_traps$has_height == 1, 19, 1) , 
       cex = 1.5, 
       col =  alpha( sex_col , 0.6 )  )
  for (i in 1:250) {
    lines(x = year_seq ,  
          y = dat_traps$min_height + post$kappa[i,1] * ( 1 - exp(-post$lambda[i,1] * year_seq)) ,  
          type = "l", 
          col = col.alpha( boycol, alpha = 0.1))}
  for (i in 1:250) {
    lines(x = year_seq ,  
          y = dat_traps$min_height + post$kappa[i,2] * ( 1 - exp(-post$lambda[i,2] * year_seq)) ,  
          type = "l", 
          col = col.alpha( girlcol, alpha = 0.1))}
  
  #plot of returns
  phi <-  mean(post_t$iota) +
    mean(post_t$gamma) * log(1-exp(- mean(post_t$beta) * seq_trait/4  )) +
    mean(post$eta_h)* mean(post$height_merged) 
  psi <-  mean(post_t$xi) * mean(log(dat_shells$duration)) 
  samp_data <- rbern(length(seq_trait/4),  
                     1 - exp (-mean(post_t$alpha) * exp(phi) * exp(psi)))
  
  actor <- d_trapppl$index_id[which(d_trapppl$anonymeID == 13212)]
  plot(jitter(seq_trait) * mean(d_trapppl$age), samp_data, 
       xlim = c(0,age_plot), ylim = c(0, 1), 
       xlab = "Age", ylab = "p capture",
       pch = 16, col = col.alpha("lawngreen", 0.2))
  #with ideal conditions (best actor, shortest time)
  # for(i in 1:150){
  #   phi <-  post_t$iota[i,actor] + 
  #     post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait)) +
  #     post$eta_h[i]*apply(post$height_merged, 1, mean)[i] 
  #   psi <-  post_t$xi[i] * log(min(dat_traps$duration))  
  #   p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
  #   lines( seq_trait * mean(d_trapppl$age),  p, 
  #          col = col.alpha(trapcol, 0.2), lwd = 1)
  # }
  #with mean actor and time
  for(i in 1:150){
    phi <-  apply(post_t$iota,1,mean )[i] +
      post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait  )) +
      post$eta_h[i]*apply(post$height_merged, 1, mean)[i] 
    
    psi <-  post_t$xi[i] * mean(log(dat_traps$duration)) 
    p <- 1 - exp (-post_t$alpha[i] * exp(phi) * exp(psi) )
    lines( seq_trait * mean(d_trapppl$age),  p, 
           col = col.alpha(trapcol, 0.2), lwd = 1)
  }
  points(jitter(d_trapppl$age[d_traps$index_id], amount = 0.5), 
         jitter(d_traps$success, amount = 0.02), 
         pch = 16, cex = ifelse(d_traps$success == 1, 0.8, 0.6), 
         col = col.alpha(othercol, ifelse(d_traps$success == 1, 0.4, 0.1)))
  


###############################################
#GRIP ESTIMATION - draft
###############################################
  post <- extract.samples(m_shells_grip)
  
  sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
  presence_col <- ifelse(dat_shells$has_grip == "1", "deepskyblue4", "darkorange3")
  year_seq <- seq(0,4,0.2)
  
  plot(x = dat_shells$age, 
       y = apply(post$grip_merged, 2, mean), 
       xlim = c(0,4),
       xlab = "Age" , 
       ylab = "grip",
       cex.lab=1.8 , 
       cex.axis=1.8 ,
       pch = ifelse(dat_shells$has_grip == 1, 19, 1) , 
       cex = 1.5, 
       col =  alpha( sex_col , 0.6 )  )
  for (i in 1:250) {
    lines(x = year_seq ,  
          y = post$min_grip[i] + post$epsilon[i,1] * ( 1 - exp(-post$upsilon[i,1] * year_seq)) ,  
          type = "l", 
          col = col.alpha( boycol, alpha = 0.1))}
  for (i in 1:250) {
    lines(x = year_seq ,  
          y = post$min_grip[i] + post$epsilon[i,2] * ( 1 - exp(-post$upsilon[i,2] * year_seq)) ,  
          type = "l", 
          col = col.alpha( girlcol, alpha = 0.1))}
  
  #plot of returns
  phi <-  mean(post$iota) +
    mean(post$gamma) * log(1-exp(- mean(post$beta) * seq_trait  )) +
    mean(post$theta_g)* mean(log(post$grip_merged)) 
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
      post$theta_g[i]* log(apply(post$grip_merged, 1, mean)[i]) 
    
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