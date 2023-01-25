
source("project_support.R")

processed_data <- list.load("data/processed_data.RData")

mean_age_shells <- mean(processed_data$shell_ppl$age)
mean_age_traps <- mean(processed_data$trap_ppl$age)

##############################
#fit to data - age only
##############################
#make data lists
dat_shells <- make_list_data_age(processed_data, foraging_type = "shells")
dat_traps <- make_list_data_age(processed_data, foraging_type = "traps")

#plot fit to data with age only
# post_s <- extract.samples(m_shell_age)
# post_t <- extract.samples(m_trap_age)

load(file = "model_fits/post_s_age.rda")
load(file = "model_fits/post_t_age.rda")
#loadRDS


png("plots/age_only.png", height = 16, width = 16, units = "cm", res = 500, type="cairo")
  par(mfrow = c(2,2),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
  #shells
  phi <-  mean(post_s$iota) +
          mean(post_s$gamma) * log(1-exp(- mean(post_s$beta) * seq_trait  )) 
  psi <-  mean(post_s$xi) * (mean(log(dat_shells$duration))) + 
          mean(post_s$tau)* mean(dat_shells$tide) 
  samp_data <- rlnorm(length(seq_trait),  
                      mean(log(post_s$alpha)) + phi + psi, 
                      mean(post_s$sigma))
  
  plot(jitter(seq_trait) * mean_age_shells, samp_data, 
       xlim = c(0,age_plot), ylim = c(0, max(dat_shells$returns)+1), 
       xlab = "Age", ylab = "kg shellfish",
       pch = 16, col = col.alpha("orange", 0.2))
  for(i in 1:150){
    phi <-  apply(post_s$iota,1,mean )[i] +
            post_s$gamma[i] * log(1-exp(- post_s$beta[i] * seq_trait  )) 
    psi <-  post_s$xi[i] * (mean(log(dat_shells$duration))) +  
            post_s$tau[i]* mean(dat_shells$tide)
    R <- exp (  log(post_s$alpha)[i] + phi + psi + 
                  (post_s$sigma[i]^2 /2))
    lines( seq_trait * mean_age_shells,  R, 
           col = col.alpha(shellcol, 0.2), lwd = 1)
  }
  points(jitter(dat_shells$age[dat_shells$ID_i] * mean_age_shells, amount = 0.4), 
         as.numeric(dat_shells$returns),
         pch = 16, cex = 0.8, col = col.alpha(othercol, 0.4))
  text(1, max(dat_shells$returns)+0.5, "A")

  #traps 
  phi <-  mean(post_t$iota) +
          mean(post_t$gamma) * log(1-exp(- mean(post_t$beta) * seq_trait  )) 
  psi <-  mean(post_t$xi) * (mean(log (dat_traps$duration))) 
  p <- 1 - exp ( - mean(post_t$alpha) * exp(phi) * exp(psi))
  samp_data <- rbern(length(seq_trait),  p)
  
  #NB making plot with 13212 only (one of best hunters) to make plot with optimal situation to raise curve off zero
  plot(jitter(seq_trait) * mean_age_traps, samp_data, 
       xlab = "Age", ylab = "p capture",
       xlim = c(0,age_plot), ylim = c(0, 1), 
       pch = 16, col = col.alpha("lawngreen", 0.3))
  #with ideal conditions (best actor, shortest time)
  for(i in 1:150){
    phi <-  post_t$iota[i,dat_traps$best_guy] + 
            post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait)) 
    psi <-  post_t$xi[i] * log(min(dat_traps$duration))  
    p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
    lines( seq_trait * mean_age_traps,  p, 
           col = col.alpha(trapcol, 0.2), lwd = 1)
  }
  #with average actor and average time 
  # for(i in 1:150){
  #   phi <-  apply(post_t$iota,1,mean )[i] + 
  #           post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait)) 
  #   psi <-  post_t$xi[i] * mean(log(dat_traps$duration))  
  #   p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
  #   lines( seq_trait * mean(d_trapppl$age),  p, 
  #          col = col.alpha(trapcol, 0.2), lwd = 1)
  # }
  points(jitter(dat_traps$age[dat_traps$ID_i] * mean_age_traps, amount = 0.5), 
         jitter(dat_traps$success, amount = 0.02), 
    pch = 16, cex = ifelse(dat_traps$success == 1, 0.8, 0.6), 
    col = col.alpha(othercol, ifelse(dat_traps$success == 1, 0.4, 0.1)))
  text(1, 0.95, "B")
  
#######################################
#age variation
#######################################
  plot(NULL, xlim = c(0,age_plot), ylim = c(0,1), 
       xlab = "Age", ylab = "proportion max foraging")#, main = "Age only"
  phi <- exp(median(post_s$gamma) * log(1-exp(-median(post_s$beta) * seq_trait  ))) 
  lines( seq_trait * mean_age_shells,  phi, col = col.alpha(shellcol, 1), lwd = 2)
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_s$gamma * log(1-exp(-post_s$beta * x ))), 0.95) )
  shade(mu_phi, seq_trait* mean_age_shells, col = col.alpha(shellcol, 0.1))
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_s$gamma * log(1-exp(-post_s$beta * x ))), 0.5) )
  shade(mu_phi, seq_trait* mean_age_shells, col = col.alpha(shellcol, 0.15))
  mu_phi <-   sapply ( seq_trait , function (x) PI ((1-exp(-post_s$beta * x )) ^ post_s$gamma, 0.3) )
  shade(mu_phi, seq_trait* mean_age_shells, col = col.alpha(shellcol, 0.15))
  for(i in 1:30){
    phi <- exp(post_s$gamma[i] * log(1-exp(-post_s$beta[i] * seq_trait )))
    lines( seq_trait * mean_age_shells,  phi, col = col.alpha(shellcol, 0.3))
  }
  text(1, 0.95, "C")
  
  
  plot(NULL, xlim = c(0,age_plot), ylim = c(0,1), 
       xlab = "Age", ylab = "proportion max foraging")#, main = "Age only"
  phi <-  exp(median(post_t$gamma) * log(1-exp(-median(post_t$beta) * seq_trait  ))) 
  lines( seq_trait * mean_age_traps,  phi, col = col.alpha(trapcol, 1), lwd = 2)
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_t$gamma * log(1-exp(-post_t$beta * x ))), 0.95) )
  shade(mu_phi, seq_trait* mean_age_traps, col = col.alpha(trapcol, 0.1))
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_t$gamma * log(1-exp(-post_t$beta * x ))), 0.5) )
  shade(mu_phi, seq_trait* mean_age_traps, col = col.alpha(trapcol, 0.15))
  mu_phi <-   sapply ( seq_trait , function (x) PI (exp(post_t$gamma * log(1-exp(-post_t$beta * x ))), 0.3) )
  shade(mu_phi, seq_trait* mean_age_traps, col = col.alpha(trapcol, 0.15))
  for(i in 1:30){
    phi <-  exp(post_t$gamma[i] * log(1-exp(-post_t$beta[i] * seq_trait )))
    lines( seq_trait * mean_age_traps,  phi, col = col.alpha(trapcol, 0.3))
  }
  text(1, 0.95, "D")
dev.off()

  
#####################################
#effect of knowledge and other stuff
#####################################
dat_shells <- make_list_data_all(processed_data, foraging_type = "shells")
dat_traps <- make_list_data_all(processed_data, foraging_type = "traps")

# post_s <- extract.samples(m_shells_all)
# post_t <- extract.samples(m_traps_all)

rm(post_s)
rm(post_t)
load(file = "model_fits/post_s.rda")
load(file = "model_fits/post_t.rda")


###########
#check max and min effect vs sds
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

diffs_shells <- diffs %>% 
  mutate( variable = fct_reorder(.f = variable, .x = kg_shellfish, .fun = mean))


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

diffs_traps <- diffs %>% 
  mutate( variable = fct_reorder(.f = variable, .x = p_success, .fun = mean))

out_shells_max <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(fill = col.alpha(shellcol, 0.6),
              colour = shellcol) +
  labs(y = "kg shellfish difference", x = "")+
  ylim(-max(dat_shells$returns), max(dat_shells$returns))+
  theme_classic()

out_traps_max <- ggplot(diffs_traps, aes(x = variable, y = p_success)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(fill = col.alpha(trapcol, 0.6),
              colour = trapcol) +
  labs(y = "p_success difference", x = "")+
  ylim(-1, 1)+
  theme_classic()

v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi)),
                                        rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi))),
                          p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
                          foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))


phi_trait_max <- ggplot(v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type , color = foraging_type) )+ 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(trim=FALSE)+
  scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c( col.alpha(shellcol, 0.6), col.alpha(trapcol, 0.6 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c(shellcol, trapcol),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  labs(x = "", y = "phi difference")+
  ylim(-10, 10)+
  theme_classic()+
  theme(legend.position="top")



png("plots/alldiffs_combined.png", height = 14, width = 16, units = "cm", res = 500, type="cairo")
plots <- align_plots(phi_trait_max, out_shells_max, align = 'v', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], out_traps_max, labels = c('B', 'C'), label_size = 12)

# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
dev.off()

#CALCULATE DIFFERENCE BETWEEN ONE STANDARD DEVIATION OF SIMULATED DATA
diffs_out <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 9))
colnames(diffs_out) <- c("sg","sh","sk","sd","st","tg","th","tk","td")
diffs_phi <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 6))
colnames(diffs_phi) <- c("sg","sh","sk","tg","th","tk")

sds <- c( sd(log(dat_shells$grip), na.rm = TRUE),
          sd(log(dat_shells$height), na.rm = TRUE), 
          sd(apply(post_s$knowledge, 2, mean)), #sd(log(dat_shells$knowledge))
          sd(log(dat_shells$duration)),
          sd(dat_shells$tide))

for (i in 1:5) {
  phi_min <-  apply(post_s$iota,1,mean )  +
    post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
    post_s$theta_g*(mean(log(dat_shells$grip), na.rm = TRUE) - ifelse(i == 1, sds[i], 0)) +
    post_s$eta_h* (mean(log(dat_shells$height), na.rm = TRUE) - ifelse(i == 2, sds[i], 0)) +
    #post_s$zeta_k* (mean(log(dat_shells$knowledge), na.rm = TRUE) - ifelse(i == 1, sds[i], 0))
    post_s$zeta_k*(mean(apply(post_s$knowledge, 2, mean)) - ifelse(i == 3, sds[i], 0))
  psi_min <-  post_s$xi * (mean(log(dat_shells$duration)) - ifelse(i == 4, sds[i], 0))+
    post_s$tau* (mean(dat_shells$tide) - ifelse(i == 5, sds[i], 0))
  min_k <- exp (log(post_s$alpha) + phi_min + psi_min +
                  (post_s$sigma^2 /2))
  phi_max <-  apply(post_s$iota,1,mean )  +
    post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
    post_s$theta_g* (mean(log(dat_shells$grip), na.rm = TRUE) + ifelse(i == 1, sds[i], 0)) +
    post_s$eta_h* (mean(log(dat_shells$height), na.rm = TRUE) + ifelse(i == 2, sds[i], 0)) +
    #post_s$zeta_k* (mean(log(dat_shells$knowledge), na.rm = TRUE) + ifelse(i == 1, sds[i], 0))
    post_s$zeta_k* (mean(apply(post_s$knowledge, 2, mean)) + ifelse(i == 3, sds[i], 0))
  psi_max <-  post_s$xi * (mean(log(dat_shells$duration)) + ifelse(i == 4, sds[i], 0))+
    post_s$tau* (mean(dat_shells$tide) + ifelse(i == 5, sds[i], 0))
  max_k <- exp (log(post_s$alpha) + phi_max + psi_max +
                  (post_s$sigma^2 /2))
  diffs_out[i] <- max_k - min_k
  if( i %in% 1:3) diffs_phi[i] <- phi_max - phi_min
}

diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out)),
                                  rep("tide", nrow(diffs_out))),
                    kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
)

diffs_shells <- diffs %>% 
  mutate( variable = fct_reorder(.f = variable, .x = kg_shellfish, .fun = mean))


sds <- c( sd(log(dat_traps$grip), na.rm = TRUE),
          sd(log(dat_traps$height), na.rm = TRUE), 
          sd(apply(post_t$knowledge, 2, mean)), #sd(log(dat_traps$knowledge))
          sd(log(dat_traps$duration)))

for (i in 1:4) {
  phi_min <-  apply(post_t$iota,1,mean )  +
    post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
    post_t$theta_g*(mean(log(dat_traps$grip), na.rm = TRUE) - ifelse(i == 1, sds[i], 0)) +
    post_t$eta_h* (mean(log(dat_traps$height), na.rm = TRUE) - ifelse(i == 2, sds[i], 0)) +
    #post_s$zeta_k* (mean(log(dat_traps$knowledge), na.rm = TRUE) - ifelse(i == 3, sds[i], 0))
    post_t$zeta_k*(mean(apply(post_t$knowledge, 2, mean)) - ifelse(i == 3, sds[i], 0))
  psi_min <-  post_t$xi * (mean(log(dat_traps$duration)) - ifelse(i == 4, sds[i], 0))
  min_k <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi_min))
  phi_max <-  apply(post_t$iota,1,mean )  +
    post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
    post_t$theta_g* (mean(log(dat_traps$grip), na.rm = TRUE) + ifelse(i == 1, sds[i], 0)) +
    post_t$eta_h* (mean(log(dat_traps$height), na.rm = TRUE) + ifelse(i == 2, sds[i], 0)) +
    #post_s$zeta_k* (mean(log(dat_traps$knowledge), na.rm = TRUE) ifelse(i == 1, sds[i], 0)) +
    post_t$zeta_k* (mean(apply(post_t$knowledge, 2, mean)) + ifelse(i == 3, sds[i], 0)) 
  psi_max <-  post_t$xi * (mean(log(dat_traps$duration)) + ifelse(i == 4, sds[i], 0))
  max_k <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi_max))
  diffs_out[i + 5] <- max_k - min_k
  if( i %in% 1:3) diffs_phi[i + 3] <- phi_max - phi_min
}

diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                  rep("grip", nrow(diffs_out)),
                                  rep("height", nrow(diffs_out)),
                                  rep("duration", nrow(diffs_out))),
                    p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
)

diffs_traps <- diffs %>% 
  mutate( variable = fct_reorder(.f = variable, .x = p_success, .fun = mean))


diffs_traps <- diffs %>% 
  mutate( variable = fct_reorder(.f = variable, .x = p_success, .fun = mean))

out_shells_sd <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(fill = col.alpha(shellcol, 0.6),
              colour = shellcol) +
  labs(y = "kg shellfish difference", x = "")+
  ylim(-4, 6)+
  theme_classic()

out_traps_sd <- ggplot(diffs_traps, aes(x = variable, y = p_success)) + 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(fill = col.alpha(trapcol, 0.6),
              colour = trapcol) +
  labs(y = "p_success difference", x = "")+
  ylim(-1, 1)+
  theme_classic()

v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi)),
                                        rep("knowledge", nrow(diffs_phi)),
                                        rep("grip", nrow(diffs_phi)),
                                        rep("height", nrow(diffs_phi))),
                          p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
                          foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))


phi_trait_sd <- ggplot(v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type , color = foraging_type) )+ 
  geom_hline(yintercept = 0, color = "grey90") +
  geom_violin(trim=FALSE)+
  scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c( col.alpha(shellcol, 0.6), col.alpha(trapcol, 0.6 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c(shellcol, trapcol),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
  labs(x = "", y = "phi difference")+
  ylim(-10, 10)+
  theme_classic()+
  theme(legend.position="top")


png("plots/alldiffs_combined_sd.png", height = 14, width = 16, units = "cm", res = 500, type="cairo")
plots <- align_plots(phi_trait_sd, out_shells_sd, align = 'v', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], out_traps_sd, labels = c('B', 'C'), label_size = 12)

# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
dev.off()

