library(rethinking)
library(rlist)
library(tidyverse)
library(ggridges)
real_data <- list.load("2_data_preparation/processed_data.RData")
seq_trait <- seq(0,3,0.001)

#########################
#DESCRIPTIVE STATISTICS 
#for M&M - non reproducible
#########################

source("../data/corrections.R")

d_trips        <- read.csv("../data/trips_sheet.csv")
d_activity     <- read.csv("../data/activity_sheet.csv")
d_returns      <- read.csv("../data/returns_sheet.csv")
d_participants <- read.csv("../data/participants_sheet.csv")

d_trips        <- correct_trips(d_trips)
d_activity     <- correct_activity(d_activity)
d_returns      <- correct_returns(d_returns)
d_participants <- correct_participants(d_participants)

d_anthropometrics <- read.csv("../data/anthropometrics.csv")

census <- read.csv("../data/census_BK_2020.csv")

d_knowledge <- list.load("../data/processed_data.RData")

#percentage of hunting trips
sum(str_detect(d_trips$trip_type, "fyuka"))/nrow(d_trips)

#n households
hhs <- unique(census$household_2019)
hhs <- hhs[- which (hhs == "BK144, NY201" |
                    hhs == "BK109, BK106" |
                    hhs == "BK134, BK161" |
                    hhs == "BK121, BK142" |
                    hhs == "inferred" |
                    hhs == "out_of_village" |
                    hhs == "BK144, BK172, BK176" |
                    hhs == "BK158, BK159" |
                    hhs == "dead" |
                    hhs == "REMOVED" |
                    hhs == "" |
                    is.na(hhs)
)]
length(hhs)

#people in village
nrow(census) - sum(census$household_2019 %in% c("inferred", "out_of_village", "dead", "REMOVED", ""))

#sample sizes
d_knowledge$N #people in the knowledge sample
sum(!is.na(d_anthropometrics$height)) #people for whom we have height
sum(!is.na(d_anthropometrics$weight)) #people for whom we have weight
sum(!is.na(d_anthropometrics$grip)) #people for whom we have grip strength

#######
#descriptive shellfish foraging
shells <- d_returns %>% 
  filter( trip_type %in% c("shell_collect", "crab_collect") ) %>% 
  select(trip, who)
#remove trips where weight has not been measured
shells <- shells %>%
  filter( !is.na(weight_g))
length(unique(shells$trip)) #number of shellfish foraging trips 
nrow(shells) #number of foraging trips/person 
length(unique(shells$who))

#descriptive sample
shell_ppl <- data.frame(ID = unique(shells$who),
                        age = rep(NA, length( unique(shells$who))),
                        sex = rep(NA, length( unique(shells$who)))
)
for(i in 1:nrow(shell_ppl)){
  shell_ppl$age[i] <- census$age[ which (census$ID == shell_ppl$ID[i])]
}
for(i in 1:nrow(shell_ppl)){
  shell_ppl$sex[i] <- census$sex[ which (census$ID == shell_ppl$ID[i])]
}
#proportion females
sum(shell_ppl$sex == "f")/nrow(shell_ppl)
#age range
min(shell_ppl$age, na.rm = T)
max(shell_ppl$age, na.rm = T)
mean(shell_ppl$age, na.rm = T)

#######
trap_trips <- d_trips$trip[ which (str_detect(d_trips$trip_type, "fyuka"))]
#n people who participated to trap trips
length(unique(d_participants$ID[which(d_participants$trip %in% trap_trips)]))
#n trips
length(trap_trips)

#other measures in wrangle data

#Traps
people_traps <- unique( na.omit( all_trap_data$ID_who_built ))
trap_ppl <- data.frame(ID = people_traps,
                       age = rep(NA, length( people_traps )),
                       sex = rep(NA, length( people_traps ))
)
for(i in 1:nrow(trap_ppl)){
  trap_ppl$age[i] <- census$age[ which (census$ID == trap_ppl$ID[i])]
}

for(i in 1:nrow(trap_ppl)){
  trap_ppl$sex[i] <- census$sex[ which (census$ID == trap_ppl$ID[i])]
}



#########################
#AGE ONLY
#########################
#prepare data age only
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

post_s <- extract.samples(m_shell_age)
post_t <- extract.samples(m_trap_age)

#generated data

#proportion max foraging by age 10
p_max_for_10 <- matrix(NA, 
                       nrow = 3, ncol = 2, 
                       dimnames = list(c("median", "5%PI", "94%PI"), 
                                       c("shells","traps")))
p_max_for_10[1,1] <- median((1-exp(-post_s$beta * 10/mean(d_shellppl$age)  )) ^ median(post_s$gamma))
p_max_for_10[2:3,1] <- PI(((1-exp(-post_s$beta * 10/mean(d_shellppl$age)  )) ^ post_s$gamma))

p_max_for_10[1,2] <- (1-exp(-median(post_t$beta) * 10/mean(d_trapppl$age)  )) ^ median(post_t$gamma)
p_max_for_10[2:3,2] <- PI((1-exp(-post_t$beta * 10/mean(d_trapppl$age)  )) ^ PI(post_t$gamma))
#WHY MEDIAN IS NOT BETWEEN PI????




######################################
# TIDE LEVELS
######################################
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells

for(i in 1:nrow(d_shells)){
  d_shells$age[i] <- d_shellppl$age[which (d_shellppl$anonymeID == d_shells$anonymeID[i])]
}

for(i in 1:nrow(d_shells)){
  d_shells$sex[i] <- d_shellppl$sex[which (d_shellppl$anonymeID == d_shells$anonymeID[i])]
}
d_tides <- d_shells [-which(d_shells$age >=15 & d_shells$sex == "m"),]
d_tides %>% filter(age >=20 ) %>% summarize( Mean = mean(tide_avg_depth))
d_tides %>% filter(age <=19 ) %>% summarize( Mean = mean(tide_avg_depth))

d_tides %>% filter(tide_avg_depth >0 ) %>% summarize( Mean = mean(age))
d_tides %>% filter(tide_avg_depth <=0 ) %>% summarize( Mean = mean(age))


##########
#SAVE DATA
##########
generated_data <- list(p_max_for_10 = p_max_for_10)

list.save(generated_data, "4_outcomes/generated_data.Rdata")
