# library(rethinking)
# library(rlist)
# library(tidyverse)
# library(ggridges)
# library(cowplot)
# real_data <- list.load("2_data_preparation/processed_data.RData")
# trapcol <- "#52b502" #"#5ca81e"
# shellcol <-  "#ea5914"
# othercol <- "grey30"#"#1482ac"
# seq_trait <- seq(0,3,0.005)
# age_plot <- 40


##############################
#fit to data - age only
##############################
#make data lists
dat_shells <- make_list_data_age(foraging_type = "shells")
dat_traps <- make_list_data_age(foraging_type = "traps")

#plot fit to data with age only
post_s <- extract.samples(m_shells_age)
post_t <- extract.samples(m_traps_age)


png("../plots/age_only.png", height = 16, width = 16, units = "cm", res = 500, type="cairo")
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
dat_shells <- make_list_data_all(foraging_type = "shells")
dat_traps <- make_list_data_all(foraging_type = "traps")

post_s <- extract.samples(m_shells_all)
post_t <- extract.samples(m_traps_all)




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



png("../plots/alldiffs_combined.png", height = 14, width = 16, units = "cm", res = 500, type="cairo")
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


png("../plots/alldiffs_combined_sd.png", height = 14, width = 16, units = "cm", res = 500, type="cairo")
plots <- align_plots(phi_trait_sd, out_shells_sd, align = 'v', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], out_traps_sd, labels = c('B', 'C'), label_size = 12)

# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
dev.off()



# 
# 
# #CALCULATE DIFFERENCE BETWEEN MINIMUM AND MAXIMUM ON OUTCOME SCALE
# diffs_out <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 9))
# colnames(diffs_out) <- c("sk","sh","sg","sd","st","tk","th","tg","td")
# diffs_phi <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 6))
# colnames(diffs_phi) <- c("sk","tk","sh","th","sg","tg")
# #SHELLS
# phi_min <-  apply(post_s$iota,1,mean )  +
#         post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#         post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#         post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#         post_s$zeta_k* min(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#         post_s$tau* mean(dat_shells$tide)
# min_k <- exp (log(post_s$alpha) + phi_min + psi +
#               (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#         post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#         post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#         post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#         post_s$zeta_k* max(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#         post_s$tau* mean(dat_shells$tide)
# max_k <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sk <- max_k - min_k
# diffs_phi$sk <- phi_max - phi_min
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* min(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# min_h <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* max(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# max_h <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sh <- max_h - min_h
# diffs_phi$sh <- phi_max - phi_min
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* min(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# min_g <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* max(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# phi_max <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# max_g <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sg <- max_g - min_g
# diffs_phi$sg <- phi_max - phi_min
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * min(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# min_d <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * max(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# max_d <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sd <- max_d - min_d
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* min(dat_shells$tide)
# min_t <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* max(dat_shells$tide)
# max_t <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$st <- min_t - max_t
# 
# diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
#                                   rep("grip", nrow(diffs_out)),
#                                   rep("height", nrow(diffs_out)),
#                                   rep("duration", nrow(diffs_out)),
#                                   rep("tide", nrow(diffs_out))),
#                     kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
# )
# 
# diffs_shells <- diffs %>%
#   mutate( variable = fct_reorder(.f = variable, .x = kg_shellfish, .fun = mean))
# # ggplot(diffs_shells, aes(x = kg_shellfish, y = variable)) +
# #   geom_vline(xintercept = 0, colour = col.alpha("grey40", 0.3))+
# #   geom_density_ridges(fill = col.alpha("cornflowerblue", 0.3),
# #                       colour = "cornflowerblue") +
# #   labs(x = "kg shellfish difference", y = "")+
# #   xlim(-0.5, 6)+
# #   theme_classic()
# 
# out_shells <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) +
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(fill = col.alpha(shellcol, 0.6),
#               colour = shellcol) +
#   labs(y = "kg shellfish difference", x = "")+
#   ylim(-4, 6)+
#   theme_classic()
# 
# # png("../plots/diffout_shell_violin.png", height = 10, width = 10, units = "cm", res = 500)
# # out_shells
# # dev.off()
# #
# #
# # png("../plots/diffout_shell_caterpillar.png", height = 10, width = 12, units = "cm", res = 500)
# # plot(NULL, xlim = c(1,5), ylim = c(-3,6),
# #      xlab = "", ylab = "kg shellfish difference",xaxt="n")
# # points(1:5, apply(diffs_out[,1:5], 2, mean), pch = 16, col = shellcol, size = 4 )
# # for(i in 1:5) lines(rep(i, 2), PI(diffs_out[,i]), col = shellcol, lwd = 4)
# # axis(1, at= c(1:5), labels= c("knowledge", "height", "grip", "duration", "tide"))
# # dev.off()
# 
# 
# 
# 
# #TRAPS
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* min(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration))
# min_k <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* max(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration))
# max_k <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
# diffs_out$tk <- (max_k - min_k)[1:length(post_s$alpha)]
# diffs_phi$tk <- (phi_max - phi_min)[1:length(post_s$alpha)]
# 
# 
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* min(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration))
# min_h <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* max(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration))
# max_h <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
# diffs_out$th <- (max_h - min_h)[1:length(post_s$alpha)]
# diffs_phi$th <- (phi_max - phi_min)[1:length(post_s$alpha)]
# 
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* min(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration))
# min_g <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* max(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration))
# max_g <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
# diffs_out$tg <- (max_g - min_g)[1:length(post_s$alpha)]
# diffs_phi$tg <- (phi_max - phi_min)[1:length(post_s$alpha)]
# 
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * min(log(dat_traps$duration))
# min_d <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 1.23  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * max(log(dat_traps$duration))
# max_d <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
# diffs_out$td <- (max_d - min_d )[1:length(post_s$alpha)]
# 
# 
# diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
#                                   rep("grip", nrow(diffs_out)),
#                                   rep("height", nrow(diffs_out)),
#                                   rep("duration", nrow(diffs_out))),
#                     p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
# )
# 
# diffs_traps <- diffs %>%
#   mutate( variable = fct_reorder(.f = variable, .x = p_success, .fun = mean))
# 
# out_traps <- ggplot(diffs_traps, aes(x = variable, y = p_success)) +
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(fill = col.alpha(trapcol, 0.6),
#               colour = trapcol) +
#   labs(y = "p_success difference", x = "")+
#   ylim(-1, 1)+
#   theme_classic()
# 
# # png("../plots/diffout_trap_violin.png", height = 10, width = 10, units = "cm", res = 500)
# # out_traps
# # dev.off()
# #
# #
# # png("../plots/diffout_trap_caterpillar.png", height = 10, width = 10, units = "cm", res = 500)
# # plot(NULL, xlim = c(1,4), ylim = c(-2,1),
# #      xlab = "", ylab = "p_success difference",xaxt="n")
# # points(1:4, apply(diffs_out[,6:9], 2, mean), pch = 16, col = trapcol, size = 4 )
# # for(i in 1:4) lines(rep(i, 2), PI(diffs_out[,i+5]), col = trapcol, lwd = 4)
# # axis(1, at= c(1:4), labels= c("knowledge", "height", "grip", "duration"))
# # dev.off()
# # ggplot(diffs, aes(x = p_success, y = variable)) +
# #   geom_vline(xintercept = 0, colour = col.alpha("grey40", 0.3))+
# #   geom_density_ridges(fill = col.alpha("cornflowerblue", 0.3),
# #                       colour = "cornflowerblue") +
# #   labs(x = "p_success difference", y = "")+
# #   xlim(-0.5, 1)+
# #   theme_classic()
# #
# #
# # png("../plots/phi_trait_caterpillar.png", height = 10, width = 10, units = "cm", res = 500)
# # plotcol <- rep(c(shellcol, trapcol), 3)
# # positions <- c(1.2, 1.8, 3.2, 3.8, 5.2, 5.8)
# # plot(NULL, xlim = c(1,6), ylim = c(-0.1,10),
# #      xlab = "", ylab = "diff in p success",xaxt="n")
# # points(positions, apply(diffs_phi, 2, mean), pch = 16, col = plotcol, size = 2 )
# # for(i in 1:6) lines(rep(positions[i], 2), PI(diffs_phi[,i]), col = plotcol[i], lwd = 2)
# # axis(1, at= c(1.5, 3.5, 5.5), labels= c("knowledge", "height", "grip"))
# # dev.off()
# 
# v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
#                                      rep("grip", nrow(diffs_phi)),
#                                      rep("height", nrow(diffs_phi)),
#                                      rep("knowledge", nrow(diffs_phi)),
#                                      rep("grip", nrow(diffs_phi)),
#                                      rep("height", nrow(diffs_phi))),
#                        p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
#                        foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))
# 
# 
# phi_trait <- ggplot(v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type , color = foraging_type) )+
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(trim=FALSE)+
#   scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c( col.alpha(shellcol, 0.6), col.alpha(trapcol, 0.6 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
#   scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c(shellcol, trapcol),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
#   labs(x = "", y = "phi difference")+
#   ylim(-10, 10)+
#   theme_classic()+
#   theme(legend.position="top")
# 
# # png("../plots/phi_trait_violin.png", height = 8, width = 14, units = "cm", res = 500)
# # phi_trait
# # dev.off()
# 
# 
# 
# png("../plots/alldiffs_combined.png", height = 14, width = 16, units = "cm", res = 500)
# plots <- align_plots(phi_trait, out_shells, align = 'v', axis = 'l')
# # then build the bottom row
# bottom_row <- plot_grid(plots[[2]], out_traps, labels = c('B', 'C'), label_size = 12)
# 
# # then combine with the top row for final plot
# plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
# dev.off()
# 
# 
# 
# post_t <- extract.samples(m_trap_age)
#  #traps 
#   phi <-  mean(post_t$iota) +
#           mean(post_t$gamma) * log(1-exp(- mean(post_t$beta) * seq_trait  )) 
#   psi <-  mean(post_t$xi) * (mean(log (dat_traps$duration))) 
#   lambda <- mean(post_t$alpha) * exp(phi) * exp(psi)
#   samp_data <- rpois(length(seq_trait),  lambda)
#   
#   #NB making plot with 13212 only (one of best hunters) to make plot with optimal situation to raise curve off zero
#   actor <- d_trapppl$index_id[which(d_trapppl$anonymeID == 13212)]
#   plot(jitter(seq_trait) * mean(d_trapppl$age), samp_data, 
#        xlab = "Age", ylab = "p capture",
#        xlim = c(0,age_plot), ylim = c(0, 3), 
#        pch = 16, col = col.alpha("lawngreen", 0.3))
#   #with ideal conditions (best actor, shortest time)
#   for(i in 1:150){
#     phi <-  post_t$iota[i,actor] + 
#             post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait)) 
#     psi <-  post_t$xi[i] * log(max(dat_traps$duration))  
#     lambda <-  post_t$alpha[i] * exp(phi) * exp(psi)
#     lines( seq_trait * mean(d_trapppl$age),  lambda, 
#            col = col.alpha(trapcol, 0.2), lwd = 1)
#   }
#   #with average actor and average time 
#   for(i in 1:150){
#     phi <-  apply(post_t$iota,1,mean )[i] +
#             post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait))
#     psi <-  post_t$xi[i] * mean(log(dat_traps$duration))
#     p <- 1 - exp ( - post_t$alpha[i] * exp(phi) * exp(psi))
#     lines( seq_trait * mean(d_trapppl$age),  p,
#            col = col.alpha(trapcol, 0.2), lwd = 1)
#   }
#   points(jitter(d_trapppl$age[d_traps$index_id], amount = 0.5), 
#          jitter(d_traps$success, amount = 0.02), 
#     pch = 16, cex = ifelse(d_traps$success >= 1, 0.8, 0.6), 
#     col = col.alpha(othercol, ifelse(d_traps$success >= 1, 0.4, 0.1)))
#   text(1, 0.95, "B")
# 
# post_s <- extract.samples(m_shells_all)
# post_t <- extract.samples(m_traps_all)
# 
# #CALCULATE DIFFERENCE BETWEEN MINIMUM AND MAXIMUM ON OUTCOME SCALE
# 
# #CALCULATE DIFFERENCE BETWEEN MINIMUM AND MAXIMUM ON OUTCOME SCALE
# diffs_out <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 9))
# colnames(diffs_out) <- c("sk","sh","sg","sd","st","tk","th","tg","td")
# diffs_phi <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 6))
# colnames(diffs_phi) <- c("sk","tk","sh","th","sg","tg")
# 
# 
# 
# 
# sds <- c( sd(log(dat_shells$height), na.rm = T), 
#           sd(log(dat_shells$grip), na.rm = T),
#           sd(apply(post_s$knowledge, 2, mean)), #sd(log(dat_shells$knowledge), na.rm = T)
#           sd(log(dat_shells$duration), na.rm = T),
#           sd(dat_shells$tide, na.rm = T))
# 
# for (i in 1:5) {
#   phi_min <-  apply(post_s$iota,1,mean )  +
#     post_s$gamma * log(1-exp(- post_s$beta * 1.6  )) +
#     post_s$eta_h* (mean(log(dat_shells$height), na.rm = TRUE) - ifelse(i == 1, sds[i], 0)) +
#     post_s$theta_g* (mean(log(dat_shells$grip), na.rm = TRUE) - ifelse(i == 2, sds[i], 0)) +
#     #post_s$zeta_k* (mean(log(dat_shells$knowledge), na.rm = TRUE) - ifelse(i == 1, sds[i], 0))
#     post_s$zeta_k* (mean(apply(post_s$knowledge, 2, mean)) - ifelse(i == 3, sds[i], 0))
#   psi_min <-  post_s$xi * (mean(log(dat_shells$duration)) - ifelse(i == 4, sds[i], 0))+
#     post_s$tau* (mean(dat_shells$tide) - ifelse(i == 5, sds[i], 0))
#   min_k <- exp (log(post_s$alpha) + phi_min + psi_min +
#                   (post_s$sigma^2 /2))
#   phi_max <-  apply(post_s$iota,1,mean )  +
#     post_s$gamma * log(1-exp(- post_s$beta * 1.6  )) +
#     post_s$eta_h* (mean(log(dat_shells$height), na.rm = TRUE) + ifelse(i == 1, sds[i], 0)) +
#     post_s$theta_g* (mean(log(dat_shells$grip), na.rm = TRUE) + ifelse(i == 2, sds[i], 0)) +
#     #post_s$zeta_k* (mean(log(dat_shells$knowledge), na.rm = TRUE) + ifelse(i == 1, sds[i], 0))
#     post_s$zeta_k* (mean(apply(post_s$knowledge, 2, mean)) + ifelse(i == 3, sds[i], 0))
#   psi_max <-  post_s$xi * (mean(log(dat_shells$duration)) + ifelse(i == 4, sds[i], 0))+
#     post_s$tau* (mean(dat_shells$tide) + ifelse(i == 5, sds[i], 0))
#   max_k <- exp (log(post_s$alpha) + phi_max + psi_max +
#                   (post_s$sigma^2 /2))
#   diffs_out[i] <- max_k - min_k
#   if( i %in% 1:3) diffs_phi[i] <- phi_max - phi_min
# }
# 
# diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
#                                   rep("grip", nrow(diffs_out)),
#                                   rep("height", nrow(diffs_out)),
#                                   rep("duration", nrow(diffs_out)),
#                                   rep("tide", nrow(diffs_out))),
#                     kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
# )
# 
# diffs_shells <- diffs %>% 
#   mutate( variable = fct_reorder(.f = variable, .x = kg_shellfish, .fun = mean))
# 
# 
# sds <- c( sd(log(dat_traps$height)), 
#           sd(log(dat_traps$grip)),
#           sd(apply(post_s$knowledge, 2, mean)), #sd(log(dat_shells$knowledge))
#           sd(log(dat_shells$duration)),
#           sd(dat_shells$tide))
# 
# for (i in 1:4) {
#   phi_min <-  apply(post_t$iota,1,mean )  +
#     post_t$gamma * log(1-exp(- post_t$beta * 1.6  )) +
#     post_t$eta_h* (mean(log(dat_traps$height), na.rm = TRUE) - ifelse(i == 1, sds[i], 0)) +
#     post_t$theta_g* (mean(log(dat_traps$grip), na.rm = TRUE) - ifelse(i == 2, sds[i], 0)) +
#     #post_s$zeta_k* (mean(log(dat_traps$knowledge), na.rm = TRUE) - ifelse(i == 3, sds[i], 0))
#     post_t$zeta_k* (mean(apply(post_t$knowledge, 2, mean)) - ifelse(i == 3, sds[i], 0))
#   psi <-  post_t$xi * (mean(log(dat_traps$duration)) - ifelse(i == 4, sds[i], 0))
#   min_k <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi))
#   phi_max <-  apply(post_t$iota,1,mean )  +
#     post_t$gamma * log(1-exp(- post_t$beta * 1.6  )) +
#     post_t$eta_h* (mean(log(dat_traps$height), na.rm = TRUE) + ifelse(i == 1, sds[i], 0)) +
#     post_t$theta_g* (mean(log(dat_traps$grip), na.rm = TRUE) + ifelse(i == 1, sds[i], 0)) +
#     #post_s$zeta_k* (mean(log(dat_traps$knowledge), na.rm = TRUE) ifelse(i == 1, sds[i], 0)) +
#     post_t$zeta_k* (mean(apply(post_t$knowledge, 2, mean)) + ifelse(i == 1, sds[i], 0)) 
#   psi <-  post_t$xi * mean(log(dat_traps$duration)) 
#   max_k <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi))
#   diffs_out$tk <- (max_k - min_k)[1:length(post_s$alpha)]
#   diffs_phi$tk <- (phi_max - phi_min)[1:length(post_s$alpha)]
#   diffs_out[i + 5] <- max_k - min_k
#   if( i %in% 1:3) diffs_phi[i + 5] <- phi_max - phi_min
# }
# 
# diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
#                                   rep("grip", nrow(diffs_out)),
#                                   rep("height", nrow(diffs_out)),
#                                   rep("duration", nrow(diffs_out))),
#                     p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
# )
# 
# diffs_traps <- diffs %>% 
#   mutate( variable = fct_reorder(.f = variable, .x = p_success, .fun = mean))
# 
# out_shells <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) + 
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(fill = col.alpha(shellcol, 0.6),
#               colour = shellcol) +
#   labs(y = "kg shellfish difference", x = "")+
#   ylim(-4, 6)+
#   theme_classic()
# 
# out_traps <- ggplot(diffs_traps, aes(x = variable, y = p_success)) + 
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(fill = col.alpha(trapcol, 0.6),
#               colour = trapcol) +
#   labs(y = "p_success difference", x = "")+
#   ylim(-1, 1)+
#   theme_classic()
# 
# v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
#                                         rep("grip", nrow(diffs_phi)),
#                                         rep("height", nrow(diffs_phi)),
#                                         rep("knowledge", nrow(diffs_phi)),
#                                         rep("grip", nrow(diffs_phi)),
#                                         rep("height", nrow(diffs_phi))),
#                           p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
#                           foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))
# 
# 
# phi_trait <- ggplot(v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type , color = foraging_type) )+ 
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(trim=FALSE)+
#   scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c( col.alpha(shellcol, 0.6), col.alpha(trapcol, 0.6 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
#   scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c(shellcol, trapcol),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
#   labs(x = "", y = "phi difference")+
#   ylim(-10, 10)+
#   theme_classic()+
#   theme(legend.position="top")
# 
# png("../plots/alldiffs_combined_newmodels.png", height = 14, width = 16, units = "cm", res = 500)
# plots <- align_plots(phi_trait, out_shells, align = 'v', axis = 'l')
# # then build the bottom row
# bottom_row <- plot_grid(plots[[2]], out_traps, labels = c('B', 'C'), label_size = 12)
# 
# # then combine with the top row for final plot
# plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
# dev.off()
# 
# 
# diffs_out <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 9))
# colnames(diffs_out) <- c("sk","sh","sg","sd","st","tk","th","tg","td")
# diffs_phi <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 6))
# colnames(diffs_phi) <- c("sk","tk","sh","th","sg","tg")
# #SHELLS
# phi_min <-  apply(post_s$iota,1,mean )  +
#         post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#         post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#         apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#         post_s$zeta_k* min(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#         post_s$tau* mean(dat_shells$tide)
# min_k <- exp (log(post_s$alpha) + phi_min + psi +
#               (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#         post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#         post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#         apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#         post_s$zeta_k* max(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#         post_s$tau* mean(dat_shells$tide)
# max_k <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sk <- max_k - min_k
# diffs_phi$sk <- phi_max - phi_min
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* min(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# min_h <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* max(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# max_h <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sh <- max_h - min_h
# diffs_phi$sh <- phi_max - phi_min
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * min(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# min_g <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * max(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# phi_max <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# max_g <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sg <- max_g - min_g
# diffs_phi$sg <- phi_max - phi_min
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * min(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# min_d <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * max(log(dat_shells$duration)) +
#   post_s$tau* mean(dat_shells$tide)
# max_d <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$sd <- max_d - min_d
# 
# phi_min <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* min(dat_shells$tide)
# min_t <- exp (log(post_s$alpha) + phi_min + psi +
#                 (post_s$sigma^2 /2))
# phi_max <-  apply(post_s$iota,1,mean )  +
#   post_s$gamma * log(1-exp(- post_s$beta * 20  )) +
#   post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
#   apply(post_s$theta_g, 1, mean) * mean(log(dat_shells$grip), na.rm = TRUE) +
#   post_s$zeta_k* mean(apply(post_s$knowledge, 2, mean))
# psi <-  post_s$xi * mean(log(dat_shells$duration)) +
#   post_s$tau* max(dat_shells$tide)
# max_t <- exp (log(post_s$alpha) + phi_max + psi +
#                 (post_s$sigma^2 /2))
# diffs_out$st <- min_t - max_t
# 
# diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
#                                   rep("grip", nrow(diffs_out)),
#                                   rep("height", nrow(diffs_out)),
#                                   rep("duration", nrow(diffs_out)),
#                                   rep("tide", nrow(diffs_out))),
#                     kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
# )
# 
# diffs_shells <- diffs %>% 
#   mutate( variable = fct_reorder(.f = variable, .x = kg_shellfish, .fun = mean))
# ggplot(diffs_shells, aes(x = kg_shellfish, y = variable)) + 
#   geom_vline(xintercept = 0, colour = col.alpha("grey40", 0.3))+
#   geom_density_ridges(fill = col.alpha("cornflowerblue", 0.3), 
#                       colour = "cornflowerblue") +
#   labs(x = "kg shellfish difference", y = "")+
#   xlim(-0.5, 6)+
#   theme_classic()
# 
# out_shells <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) + 
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(fill = col.alpha(shellcol, 0.6),
#               colour = shellcol) +
#   labs(y = "kg shellfish difference", x = "")+
#   ylim(-4, 6)+
#   theme_classic()
# 
#   
#   #TRAPS
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* min(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration)) 
# min_k <- post_t$alpha * exp(phi_min) * exp(psi)
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* max(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration)) 
# max_k <-  post_t$alpha * exp(phi_max) * exp(psi)
# diffs_out$tk <- (max_k - min_k)[1:length(post_s$alpha)]
# diffs_phi$tk <- (phi_max - phi_min)[1:length(post_s$alpha)]
# 
# 
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* min(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration)) 
# min_h <-  post_t$alpha * exp(phi_min) * exp(psi)
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* max(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration)) 
# max_h <- post_t$alpha * exp(phi_max) * exp(psi)
# diffs_out$th <- (max_h - min_h)[1:length(post_s$alpha)]
# diffs_phi$th <- (phi_max - phi_min)[1:length(post_s$alpha)]
# 
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* min(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration)) 
# min_g <-  post_t$alpha * exp(phi_min) * exp(psi)
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* max(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * mean(log(dat_traps$duration)) 
# max_g <-  post_t$alpha * exp(phi_max) * exp(psi)
# diffs_out$tg <- (max_g - min_g)[1:length(post_s$alpha)]
# diffs_phi$tg <- (phi_max - phi_min)[1:length(post_s$alpha)]
# 
# phi_min <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * min(log(dat_traps$duration)) 
# min_d <-  post_t$alpha * exp(phi_min) * exp(psi)
# phi_max <-  apply(post_t$iota,1,mean )  +
#   post_t$gamma * log(1-exp(- post_t$beta * 20  )) +
#   post_t$eta_h* mean(log(dat_traps$height), na.rm = TRUE) +
#   post_t$theta_g* mean(log(dat_traps$grip), na.rm = TRUE) +
#   post_t$zeta_k* mean(apply(post_t$knowledge, 2, mean))
# psi <-  post_t$xi * max(log(dat_traps$duration)) 
# max_d <- post_t$alpha * exp(phi_max) * exp(psi)
# diffs_out$td <- (max_d - min_d )[1:length(post_s$alpha)]
# 
# 
# diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
#                                   rep("grip", nrow(diffs_out)),
#                                   rep("height", nrow(diffs_out)),
#                                   rep("duration", nrow(diffs_out))),
#                     p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
# )
# 
# diffs_traps <- diffs %>% 
#   mutate( variable = fct_reorder(.f = variable, .x = p_success, .fun = mean))
# 
# out_traps <- ggplot(diffs_traps, aes(x = variable, y = p_success)) + 
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(fill = col.alpha(trapcol, 0.6),
#               colour = trapcol) +
#   labs(y = "p_success difference", x = "")+
#   ylim(-1, 1)+
#   theme_classic()
# 
# v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
#                                      rep("grip", nrow(diffs_phi)),
#                                      rep("height", nrow(diffs_phi)),
#                                      rep("knowledge", nrow(diffs_phi)),
#                                      rep("grip", nrow(diffs_phi)),
#                                      rep("height", nrow(diffs_phi))),
#                        p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
#                        foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))
# 
# 
# phi_trait <- ggplot(v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type , color = foraging_type) )+ 
#   geom_hline(yintercept = 0, color = "grey90") +
#   geom_violin(trim=FALSE)+
#   scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c( col.alpha(shellcol, 0.6), col.alpha(trapcol, 0.6 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
#   scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c(shellcol, trapcol),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
#   labs(x = "", y = "phi difference")+
#   ylim(-10, 10)+
#   theme_classic()+
#   theme(legend.position="top")
# 
# png("../plots/alldiffs_combined_newmodels.png", height = 14, width = 16, units = "cm", res = 500)
# plots <- align_plots(phi_trait, out_shells, align = 'v', axis = 'l')
# # then build the bottom row
# bottom_row <- plot_grid(plots[[2]], out_traps, labels = c('B', 'C'), label_size = 12)
# 
# # then combine with the top row for final plot
# plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
# dev.off()
