
source("project_support.R")

tide_data <- read.csv("data/tide_data.csv")
processed_data <- list.load("data/processed_data.RData")

mean_age_shells <- mean(processed_data$shell_ppl$age)
mean_age_traps <- mean(processed_data$trap_ppl$age)


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

png("plots/prior_predictive_simulation.png", height = 3, width = 5, units = "in", res = 500, type="cairo")
#plot prior predictive simulation
par(mfrow = c(1,2),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
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
#dev.off()


# png("plots/prior_age_and_trait_4.png", height = 3, width = 4, units = "in", res = 500)
# #plot prior predictive simulation
#par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(NULL, xlim = c(0,3), ylim = c (0,2), 
     xlab = "Age", ylab = "Proportion improvement",
     xaxt="n")
axis(1, at=seq(0,3,0.5),labels=seq(0,30,5))
#calculate per sample
lines (AGE, (1-exp(-1 * AGE  )) ^ 1,
       col = col.alpha("darkblue", 0.7))
text(0.5, 1.9,expression(paste(beta,"=1,",gamma,"=1",sep = "")))
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
 # lines( AGE,  phi[i,], col = col.alpha("lightblue", 0.3))
}

############################################################
#CALCULATE PRIORS FOR EFFECTS OF KNOWLEDGE ETC
############################################################
#load data
dat_shells <- make_list_data_all(processed_data, foraging_type = "shells")
dat_traps <- make_list_data_all(processed_data, foraging_type = "traps")

#load fit to extract knowledge values and overlay posterior
# post_s <- extract.samples(m_shells_all)
# post_t <- extract.samples(m_traps_all)

load(file = "model_fits/post_s.rda")
load(file = "model_fits/post_t.rda")

#define priors
n_samp <- nrow(post_s$alpha)
iota <- rnorm(n_samp, 0, 1)
alpha <- rtruncnorm(n_samp, 0, 1)
beta <- rexp(n_samp, 1)
gamma <- rexp(n_samp, 1)
zeta_k <- rnorm(n_samp, 0, 1)
eta_h <- rnorm(n_samp, 0, 1)
theta_g <- rnorm(n_samp, 0, 1)
xi <- rnorm(n_samp, 0, 1)
tau <- rnorm(n_samp, 0, 1)
sigma <- rexp(n_samp, 1)

###########
#prior
###########

diffs_out <- data.frame(matrix(nrow=length(alpha), ncol = 9))
colnames(diffs_out) <- c("sg","sh","sk","sd","st","tg","th","tk","td")
diffs_phi <- data.frame(matrix(nrow=length(alpha), ncol = 6))
colnames(diffs_phi) <- c("sg","sh","sk","tg","th","tk")

#max - min - prior
for (i in 1:5) {
  phi_min <-  iota  +
    gamma * log(1-exp(- beta * 1.23  )) +
    theta_g* ifelse(i == 1 , min(log(dat_shells$grip), na.rm = TRUE), mean(log(dat_shells$grip), na.rm = TRUE))  +
    eta_h* ifelse(i == 2 , min(log(dat_shells$height), na.rm = TRUE), mean(log(dat_shells$height), na.rm = TRUE)) +
    #zeta_k* * ifelse(i == 3, min(log(dat_shells$knowledge), na.rm = TRUE) , mean(log(dat_shells$knowledge), na.rm = TRUE))
    zeta_k* ifelse(i == 3 , min(apply(post_s$knowledge, 2, mean)) , mean(apply(post_s$knowledge, 2, mean)) )
  psi_min <-  xi * ifelse(i == 4 , min(log(dat_shells$duration), na.rm = TRUE), mean(log(dat_shells$duration), na.rm = TRUE))+
    tau*ifelse(i == 5 , min(dat_shells$tide, na.rm = TRUE), mean(dat_shells$tide, na.rm = TRUE))
  min_r <- exp (log(alpha) + phi_min + psi_min +
                  (sigma^2 /2))
  phi_max <-  iota  +
    gamma * log(1-exp(- beta * 1.23  )) +
    theta_g* ifelse(i == 1 , max(log(dat_shells$grip), na.rm = TRUE), mean(log(dat_shells$grip), na.rm = TRUE))  +
    eta_h* ifelse(i == 2 , max(log(dat_shells$height), na.rm = TRUE), mean(log(dat_shells$height), na.rm = TRUE)) +
    #zeta_k* * ifelse(i == 3, max(log(dat_shells$knowledge), na.rm = TRUE) , mean(log(dat_shells$knowledge), na.rm = TRUE))
    zeta_k* ifelse(i == 3 , max(apply(post_s$knowledge, 2, mean)) , mean(apply(post_s$knowledge, 2, mean)) )
  psi_max <-  xi * ifelse(i == 4 , max(log(dat_shells$duration), na.rm = TRUE), mean(log(dat_shells$duration), na.rm = TRUE))+
    tau*ifelse(i == 5 , max(dat_shells$tide, na.rm = TRUE), mean(dat_shells$tide, na.rm = TRUE))
  max_r <- exp (log(alpha) + phi_max + psi_max +
                  (sigma^2 /2))
  diffs_out[i] <- max_r - min_r
  if( i %in% 1:3) diffs_phi[i] <- phi_max - phi_min
}
diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out)),
                                  rep("tide", nrow(diffs_out))),
                    kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
)

diffs_shells_prior <- diffs 

for (i in 1:4) {
  phi_min <-  iota  +
    gamma * log(1-exp(- beta * 1.23  )) +
    theta_g* ifelse( i == 1, min(log(dat_traps$grip), na.rm = TRUE), mean(log(dat_traps$grip), na.rm = TRUE) )  +
    eta_h* ifelse( i == 2, min(log(dat_traps$height), na.rm = TRUE), mean(log(dat_traps$height), na.rm = TRUE) ) +
    #zeta_k* ifelse( i == 3, min(log(dat_traps$knowledge), na.rm = TRUE), mean(log(dat_traps$knowledge), na.rm = TRUE) )
    zeta_k* ifelse( i == 3, min(apply(post_t$knowledge, 2, mean)), mean(apply(post_t$knowledge, 2, mean) ) ) 
  psi <-  xi * ifelse( i == 4, min(log(dat_traps$duration), na.rm = TRUE), mean(log(dat_traps$duration), na.rm = TRUE) ) 
  min_k <- 1 - exp ( - alpha * exp(phi_min) * exp(psi))
  phi_max <-  iota  +
    gamma * log(1-exp(- beta * 1.23  )) +
    theta_g* ifelse( i == 1, max(log(dat_traps$grip), na.rm = TRUE), mean(log(dat_traps$grip), na.rm = TRUE) ) +
    eta_h * ifelse( i == 2, max(log(dat_traps$height), na.rm = TRUE), mean(log(dat_traps$height), na.rm = TRUE) ) +
    #zeta_k* ifelse( i == 3, max(log(dat_traps$knowledge), na.rm = TRUE), mean(log(dat_traps$knowledge), na.rm = TRUE) ) +
    zeta_k* ifelse( i == 3, max(apply(post_t$knowledge, 2, mean)), mean( apply(post_t$knowledge, 2, mean) ) ) 
  psi <-  xi * ifelse( i == 4, max(log(dat_traps$height), na.rm = TRUE), mean(log(dat_traps$height), na.rm = TRUE) ) 
  max_k <- 1 - exp ( - alpha * exp(phi_max) * exp(psi))
  diffs_out[i + 5] <- max_k - min_k
  if( i %in% 1:3) diffs_phi[i + 3] <- phi_max - phi_min
}

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


###########
#max - min - posterior
###########

diffs_out <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 9))
colnames(diffs_out) <- c("sg","sh","sk","sd","st","tg","th","tk","td")
diffs_phi <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 6))
colnames(diffs_phi) <- c("sg","sh","sk","tg","th","tk")

#max - min - gives same result as the long one
for (i in 1:5) {
  phi_min <-  apply(post_s$iota,1,mean )  +
    post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
    post_s$theta_g* ifelse(i == 1 , min(log(dat_shells$grip), na.rm = TRUE), mean(log(dat_shells$grip), na.rm = TRUE))  +
    post_s$eta_h* ifelse(i == 2 , min(log(dat_shells$height), na.rm = TRUE), mean(log(dat_shells$height), na.rm = TRUE)) +
    #post_s$zeta_k* * ifelse(i == 3, min(log(dat_shells$knowledge), na.rm = TRUE) , mean(log(dat_shells$knowledge), na.rm = TRUE))
    post_s$zeta_k* ifelse(i == 3 , min(apply(post_s$knowledge, 2, mean)) , mean(apply(post_s$knowledge, 2, mean)) )
  psi_min <-  post_s$xi * ifelse(i == 4 , min(log(dat_shells$duration), na.rm = TRUE), mean(log(dat_shells$duration), na.rm = TRUE))+
    post_s$tau*ifelse(i == 5 , min(dat_shells$tide, na.rm = TRUE), mean(dat_shells$tide, na.rm = TRUE))
  min_r <- exp (log(post_s$alpha) + phi_min + psi_min +
                  (post_s$sigma^2 /2))
  phi_max <-  apply(post_s$iota,1,mean )  +
    post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
    post_s$theta_g* ifelse(i == 1 , max(log(dat_shells$grip), na.rm = TRUE), mean(log(dat_shells$grip), na.rm = TRUE))  +
    post_s$eta_h* ifelse(i == 2 , max(log(dat_shells$height), na.rm = TRUE), mean(log(dat_shells$height), na.rm = TRUE)) +
    #post_s$zeta_k* * ifelse(i == 3, max(log(dat_shells$knowledge), na.rm = TRUE) , mean(log(dat_shells$knowledge), na.rm = TRUE))
    post_s$zeta_k* ifelse(i == 3 , max(apply(post_s$knowledge, 2, mean)) , mean(apply(post_s$knowledge, 2, mean)) )
  psi_max <-  post_s$xi * ifelse(i == 4 , max(log(dat_shells$duration), na.rm = TRUE), mean(log(dat_shells$duration), na.rm = TRUE))+
    post_s$tau*ifelse(i == 5 , max(dat_shells$tide, na.rm = TRUE), mean(dat_shells$tide, na.rm = TRUE))
  max_r <- exp (log(post_s$alpha) + phi_max + psi_max +
                  (post_s$sigma^2 /2))
  diffs_out[i] <- max_r - min_r
  if( i %in% 1:3) diffs_phi[i] <- phi_max - phi_min
}
diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out)),
                                  rep("tide", nrow(diffs_out))),
                    kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
)

diffs_shells <- diffs


for (i in 1:4) {
  phi_min <-  apply(post_t$iota,1,mean )  +
    post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
    post_t$theta_g* ifelse( i == 1, min(log(dat_traps$grip), na.rm = TRUE), mean(log(dat_traps$grip), na.rm = TRUE) )  +
    post_t$eta_h* ifelse( i == 2, min(log(dat_traps$height), na.rm = TRUE), mean(log(dat_traps$height), na.rm = TRUE) ) +
    #post_s$zeta_k* ifelse( i == 3, min(log(dat_traps$knowledge), na.rm = TRUE), mean(log(dat_traps$knowledge), na.rm = TRUE) )
    post_t$zeta_k* ifelse( i == 3, min(apply(post_t$knowledge, 2, mean)), mean(apply(post_t$knowledge, 2, mean) ) ) 
  psi <-  post_t$xi * ifelse( i == 4, min(log(dat_traps$duration), na.rm = TRUE), mean(log(dat_traps$duration), na.rm = TRUE) ) 
  min_k <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
  phi_max <-  apply(post_t$iota,1,mean )  +
    post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
    post_t$theta_g* ifelse( i == 1, max(log(dat_traps$grip), na.rm = TRUE), mean(log(dat_traps$grip), na.rm = TRUE) ) +
    post_t$eta_h * ifelse( i == 2, max(log(dat_traps$height), na.rm = TRUE), mean(log(dat_traps$height), na.rm = TRUE) ) +
    #post_s$zeta_k* ifelse( i == 3, max(log(dat_traps$knowledge), na.rm = TRUE), mean(log(dat_traps$knowledge), na.rm = TRUE) ) +
    post_t$zeta_k* ifelse( i == 3, max(apply(post_t$knowledge, 2, mean)), mean( apply(post_t$knowledge, 2, mean) ) ) 
  psi <-  post_t$xi * ifelse( i == 4, max(log(dat_traps$height), na.rm = TRUE), mean(log(dat_traps$height), na.rm = TRUE) ) 
  max_k <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
  diffs_out[i + 5] <- max_k - min_k
  if( i %in% 1:3) diffs_phi[i + 3] <- phi_max - phi_min
}

diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out))),
                    p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
)

diffs_traps <- diffs 


v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi)),
                                        rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi))),
                          p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
                          foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))




out_shells <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(data = diffs_shells_prior, aes(x = variable, y = kg_shellfish), fill = col.alpha("orange", 0.4),
              colour = "orange") +
  geom_violin(fill = col.alpha(shellcol, 0.4),
              colour = shellcol) +
  labs(y = "kg shellfish difference", x = "")+
  ylim(-6, 6)+
  theme_classic()


out_traps <- ggplot(diffs_traps, aes(x = variable, y = p_success)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(data = diffs_traps_prior, aes(x = variable, y = p_success), fill = col.alpha("lawngreen", 0.4),
              colour = "lawngreen") +
  geom_violin(fill = col.alpha(trapcol, 0.3),
              colour = trapcol) +
  labs(y = "p_success difference", x = "")+
  ylim(-1, 1)+
  theme_classic()


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


png("plots/alldiffs_combined_prior.png", height = 14, width = 16, units = "cm", res = 500, type="cairo")
plots <- align_plots(phi_trait_prior, out_shells, align = 'v', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], out_traps, labels = c('B', 'C'), label_size = 12)

# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
dev.off()


##################################################
#AGE ONLY - traps with average 
##################################################
# post_t <- extract.samples(m_trap_age)

load("model_fits/post_t_age.rda")

png("plots/age_only_trap_mean.png", height = 8, width = 8, units = "cm", res = 500, type="cairo")
par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)

#traps 
phi <-  mean(post_t$iota) +
  mean(post_t$gamma) * log(1-exp(- mean(post_t$beta) * seq_trait  )) 
psi <-  mean(post_t$xi) * (mean(log (dat_traps$duration))) 
p <- 1 - exp ( - mean(post_t$alpha) * exp(phi) * exp(psi))
samp_data <- rbern(length(seq_trait),  p)

plot(jitter(seq_trait) * mean(processed_data$trap_ppl$age), samp_data, 
     xlab = "Age", ylab = "p capture",
     xlim = c(0,age_plot), ylim = c(0, 1), 
     pch = 16, col = col.alpha("lawngreen", 0.3))
#with average actor and average time 
for(i in 1:150){
  phi <-  apply(post_t$iota,1,mean )[i] +
          post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait))
  psi <-  post_t$xi[i] * mean(log(dat_traps$duration))
  p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
  lines( seq_trait * mean_age_traps,  p,
         col = col.alpha(trapcol, 0.2), lwd = 1)
}
points(jitter(dat_traps$age[dat_traps$ID_i] * mean_age_traps, amount = 0.5), 
       jitter(dat_traps$success, amount = 0.02), 
       pch = 16, cex = ifelse(dat_traps$success == 1, 0.8, 0.6), 
       col = col.alpha(othercol, ifelse(dat_traps$success == 1, 0.4, 0.1)))
dev.off()

##########################################################################
#AGE ONLY -TIDE MODEL
##########################################################################
dat_shells_age <- make_list_data_age(processed_data, foraging_type = "shells") # BAB: new? what is dat_Shells_AGe?

dat_tides <- dat_shells_age [ c("M", "ID_i", "tide", "age")]

# dat_tides$age <- d_tides$age[d_tides$ID_i] # bab: ????

m_tide <- cstan(file = "models/tide_age.stan", data = dat_tides)
post_tide <- extract.samples (m_tide)

dat_tides$age <- dat_tides$age * mean_age_shells

png("plots/age_and_tide.png", height = 3, width = 4, units = "in", res = 500, type="cairo")
plot(jitter(dat_tides$age, 2) , dat_tides$tide, pch = 16, col = "#1482ac")
for ( i in 1:250) lines( seq(0,4, 0.1) * mean_age_shells, 
                         post_tide$alpha[i] + post_tide$beta[i] * seq(0,4, 0.1),
                         col = col.alpha("#1482ac", 0.3))
dev.off()
##################################################
#N ITEMS
##################################################
dat_shells <- make_list_data_age(processed_data, foraging_type = "shells")

dat_tides <- dat_shells [ c("M", "ID_i", "tide", "age")]
dat_tides$age <- d_tides$age[d_tides$ID_i]
dat_tides$age <- dat_tides$age * mean_age_shells
dat_tides$n_item_types <- processed_data$shells$n_item_types

png("plots/n_items_age.png", height = 3, width = 4, units = "in", res = 500, type="cairo")
plot(dat_tides$age[which(dat_tides$n_item_types > 0)] , 
     dat_tides$n_item_types[which(dat_tides$n_item_types > 0)],
     xlim = c(5,age_plot ), ylim = c(0, 10),
     xlab = "Age", ylab = "n types of shellfish",
     pch = 16, col = "#1482ac"  )
dev.off()


  
###################
#missing data
###################
#make plots
dat_shells <- make_list_data_all(processed_data, foraging_type = "shells")
dat_traps <- make_list_data_all(processed_data, foraging_type = "traps")

# post_s <- extract.samples(m_shells_all)
# post_t <- extract.samples(m_traps_all)

load(file = "model_fits/post_s.rda")
load(file = "model_fits/post_t.rda")

#check traits by age shells
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_knowledge == "1", shellcol, "deepskyblue4")
png("plots/missing_data_validation_shells.png", height = 5, width = 8, units = "in", res = 500, type="cairo")
par(mfrow = c(1,3),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(x = dat_shells$age * mean_age_shells, 
     y = apply(post_s$height_merged, 2, mean), 
     xlim = c(0,age_plot),
     xlab = "Age" , 
     ylab = "Height",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_height == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq  * mean_age_shells,  
        y = dat_shells$min_height + post_s$kappa[i,1] * ( 1 - exp(-post_s$lambda[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq  * mean_age_shells,  
        y = dat_shells$min_height + post_s$kappa[i,2] * ( 1 - exp(-post_s$lambda[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_shells$age  * mean_age_shells, 
     y = apply(post_s$grip_merged, 2, mean), 
     xlim = c(0,age_plot),
     xlab = "Age" , 
     ylab = "Grip",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_grip == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq* mean_age_shells,  
        y = post_s$epsilon[i,1] * ( 1 - exp(-post_s$upsilon[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq * mean_age_shells,  
        y = post_s$epsilon[i,2] * ( 1 - exp(-post_s$upsilon[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_shells$age * mean_age_shells, 
     y = apply(post_s$knowledge, 2, median), 
     xlim = c(0,age_plot),
     xlab = "Age" , 
     ylab = "Estimated knowledge",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq * mean_age_shells,  
        y = post_s$omega[i] + post_s$chi[i,1] * ( 1 - exp(-post_s$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq * mean_age_shells,  
        y = post_s$omega[i] + post_s$chi[i,2] * ( 1 - exp(-post_s$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}
dev.off()



#make plots
sex_col <- ifelse(dat_traps$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_traps$has_knowledge == "1", "deepskyblue4", "darkorange3")

png("plots/missing_data_validation_traps.png", height = 5, width = 8, units = "in", res = 500, type="cairo")
par(mfrow = c(1,3),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot(x = dat_traps$age * mean_age_traps, 
     y = apply(post_t$height_merged, 2, mean), 
     xlim = c(0,age_plot),
     xlab = "Age" , 
     ylab = "Height",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_traps$has_height == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq * mean_age_traps,  
        y = dat_traps$min_height + post_t$kappa[i,1] * ( 1 - exp(-post_t$lambda[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq * mean_age_traps,  
        y = dat_traps$min_height + post_t$kappa[i,2] * ( 1 - exp(-post_t$lambda[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_traps$age * mean_age_traps, 
     y = apply(post_t$grip_merged, 2, mean), 
     xlim = c(0,age_plot),
     xlab = "Age" , 
     ylab = "grip",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_traps$has_grip == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq * mean_age_traps,  
        y = post_t$epsilon[i,1] * ( 1 - exp(-post_t$upsilon[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq * mean_age_traps,  
        y = post_t$epsilon[i,2] * ( 1 - exp(-post_t$upsilon[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_traps$age * mean_age_traps,
     y = apply(post_t$knowledge, 2, median), 
     xlim = c(0,age_plot),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_traps$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq * mean_age_traps,
        y = post_t$omega[i] + post_t$chi[i,1] * ( 1 - exp(-post_t$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq * mean_age_traps, 
        y = post_t$omega[i] + post_t$chi[i,2] * ( 1 - exp(-post_t$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}
dev.off()

##########################
#TIDE PLOTS
##########################
png("plots/tide_distance_kernel.png", height = 10, width = 12, units = "cm", res = 500, type="cairo")
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


png("plots/tide_height_trip.png", height = 10, width = 12, units = "cm", res = 500, type="cairo")
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

png("plots/tide_height&avg.png", height = 10, width = 12, units = "cm", res = 500, type="cairo")
  ggplot(tide_data, aes(x = tide_height, y = avg_tide_depth, col = trip_length)) +
  geom_point( )+
  labs(x = "max tide height", y = "average tide height", col = "Duration (min)")+
  theme_classic()
dev.off()

png("plots/simple_data.png", height = 8, width = 16, units = "cm", res = 500, type="cairo")
par(mfrow = c(1,2),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)

plot(jitter(dat_shells$age[dat_shells$ID_i] * mean_age_shells, amount = 0.4), 
     as.numeric(dat_shells$returns), 
     xlim = c(0,age_plot), ylim = c(0, max(dat_shells$returns)+1), 
     xlab = "umri", ylab = "kg komba",
     pch = 16, col = col.alpha(shellcol, 0.8))
phi <-  mean(post_s$iota) +
  mean(post_s$gamma) * log(1-exp(- mean(post_s$beta) * seq_trait  )) 
psi <-  mean(post_s$xi) * (mean(log(dat_shells$duration))) +  
  mean(post_s$tau)* mean(dat_shells$tide)
R <- exp (  log(mean(post_s$alpha)) + phi + psi + 
              (mean(post_s$sigma)^2 /2))
lines( seq_trait * mean_age_shells,  R, 
       col = col.alpha("orange", 0.8), lwd = 3)

plot(jitter(dat_traps$age[dat_traps$ID_i] * mean_age_traps, amount = 0.5), 
     jitter(dat_traps$success, amount = 0.02),
     xlab = "umri", ylab = "p kuwinda",
     xlim = c(0,age_plot), ylim = c(0, 1), 
     pch = 16, col = col.alpha(trapcol, 0.6))
phi <-  apply(post_t$iota,1,mean )[i] + 
     post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait)) 
psi <-  post_t$xi[i] * mean(log(dat_traps$duration))  
p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
lines( seq_trait * mean_age_traps,  p, 
     col = col.alpha("lawngreen", 0.8), lwd = 3)
text(10, 0.1, paste(length(unique(processed_data$traps$trap_ID)), "fyuka"))
text(30, 0.9, paste(sum(dat_traps$success), "mawindo"))

dev.off()
