#NB THIS SCRIPT IS NOT EXECUTABLE IN ABSENCE OF NON SHARED FOLDER /data
#this script wrangles the original data, not available online, and restitutes three files with anonymized, clean data ready for analysis
#'2_Data_preparation/processed_data.RData' is the main data file for use in the analysis
#'2_data_preparation/tide_data.csv' is a small dataset used for some additional analyses of tide
#"2_data_preparation/generated_quantities.RData" contains some measures calculated from raw data that are relevant metadata and might be included in the manuscript

library(tidyverse)
library(rethinking)
library(rlist)

#########################
#LOAD FUNCTIONS AND DATA---------------------------------------------------------------------------------------------
#########################
anonyme <- function( df_column) {
  anonyme_id <- vector( length = length(df_column) )
  for(i in 1:length(df_column)){
    anonyme_id[i] <- paste( which(letters == tolower(substring(df_column[i], 1, 1))), parse_number(df_column[i]), sep = "")
  }
  return(anonyme_id)
}

source("../data/corrections.R")

d_trips        <- read.csv("../data/trips_sheet.csv")
d_activity     <- read.csv("../data/activity_sheet.csv")
d_returns      <- read.csv("../data/returns_sheet.csv")
d_participants <- read.csv("../data/participants_sheet.csv")

#correct original misspelling and other similar problems
d_trips        <- correct_trips(d_trips)
d_activity     <- correct_activity(d_activity)
d_returns      <- correct_returns(d_returns)
d_participants <- correct_participants(d_participants)

d_anthropometrics <- read.csv("../data/anthropometrics.csv")

census <- read.csv("../data/census_BK_2020.csv")

d_knowledge <- list.load("../data/processed_data.RData")


#########################
#CALCUCLATE QUANTITIES FOR MANUSCRIPT-----------------------------------------------------------------------------
##########################create list to store metadata for the manuscript
generated_quantities <- list()

#n and percentage of hunting trips
generated_quantities$p_fyuka_trips <- sum(str_detect(d_trips$trip_type, "fyuka"))/nrow(d_trips)

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
generated_quantities$n_households <- length(hhs)

#people in village
generated_quantities$n_people_bk <- nrow(census) - sum(census$household_2019 %in% c("inferred", "out_of_village", "dead", "REMOVED", ""))

#sample sizes
generated_quantities$s_size_knowledge <- d_knowledge$N #people in the knowledge sample
generated_quantities$s_size_height <- sum(!is.na(d_anthropometrics$height)) #people for whom we have height
generated_quantities$s_size_weight <- sum(!is.na(d_anthropometrics$weight)) #people for whom we have weight
generated_quantities$s_size_grip <- sum(!is.na(d_anthropometrics$grip)) #people for whom we have grip strength


###########################
#LIST OF SHELLFISH RETURNS------------------------------------------------------------------------------------------
###########################
#select trips for shell collection only
shells <- d_returns %>% 
  # option to include trip for collecting crabs that two older boys performed alone. 
  #filter( trip_type %in% c("shell_collect", "crab_collect") ) %>% 
  filter( trip_type == "shell_collect" ) %>% 
  select(entry_ID, trip, weight_g, what, who, time_start, time_end, trip_length)
#remove trips where weight has not been measured
shells <- shells %>%
  filter( !is.na(weight_g))

#merge separate datapoints for trip 679FF by person
#find number of persons and make empty dataframe
n_ppl_679FF <- length (unique (shells [shells$trip == "679FF", "who"] ))
trip_679FF <- data.frame( entry_ID = rep ("combined_entry", n_ppl_679FF),
                          trip = rep("679FF", n_ppl_679FF),
                          weight_g = rep (NA, n_ppl_679FF),
                          what = rep (NA, n_ppl_679FF),
                          who = unique (shells [shells$trip == "679FF", "who"] ),
                          time_start = rep (NA, n_ppl_679FF),
                          time_end = rep (NA, n_ppl_679FF),
                          trip_length = rep (NA, n_ppl_679FF) )
#populate dataframe
for(i in 1:nrow(trip_679FF)){
  trip_679FF$weight_g[i] <- sum( as.numeric(shells$weight_g[shells$trip == "679FF" & shells$who == trip_679FF$who[i]]))
  trip_679FF$what[i] <- paste(shells$what[shells$trip == "679FF" & shells$who == trip_679FF$who[i]], collapse = ", ")
  trip_679FF$time_start[i] <- min (shells$time_start[shells$trip == "679FF" & shells$who == trip_679FF$who[i]])
  trip_679FF$time_end[i] <- max (shells$time_end[shells$trip == "679FF" & shells$who == trip_679FF$who[i]])
  trip_679FF$trip_length[i] <- sum (shells$trip_length[shells$trip == "679FF" & shells$who == trip_679FF$who[i]])
}
#merge with main and remove other datapoints
shells <- shells %>% filter( !trip == "679FF")
shells <- rbind (trip_679FF, shells)
rm(trip_679FF)


#add date of trip
for (i in 1:nrow(shells)){
  shells$date[i] <- d_trips[ which ( d_trips$trip == shells$trip[i]), "date"]
}

##########################
#TIDE
##########################
#add tide info to shells
shells$tide_height <- NA
shells$tide_time <- NA
for (i in 1:nrow(shells)){
  shells$tide_height[i] <- as.numeric(d_trips[ which ( d_trips$trip == shells$trip[i]), "tide_height"])
  shells$tide_time[i] <- d_trips[ which ( d_trips$trip == shells$trip[i]), "tide_time"]
}
#missing tides - calculate back time of minimum and estimated tide depth
shells$tide_time[which (shells$trip == "325JD")] <- "12:52" #tide time two days earlier at 11:14, tide moves of 50 min a day on average
shells$tide_time[which (shells$trip == "009QT")] <- "10:00" #tide time following day at 10:49, tide moves of 50 min a day on average
shells$tide_height[which (shells$trip == "325JD")] <- 0.1 #tide height two days earlier at 0.15 in decreasing trend
shells$tide_height[which (shells$trip == "009QT")] <- 0.15 #tide time following day at 0.14 in slowly decreasing trend

#define format dates
shells$start_time <- paste(shells$date, " ", shells$time_start, sep = "")
shells$end_time   <- paste(shells$date, " ", shells$time_end, sep = "")
shells$tide_time  <- paste(shells$date, " ", shells$tide_time, sep = "")

shells$start_time <- as.POSIXct(shells$start_time, format="%Y-%m-%d %H:%M")
shells$end_time   <- as.POSIXct(shells$end_time, format="%Y-%m-%d %H:%M")
shells$tide_time  <- as.POSIXct(shells$tide_time, format="%Y-%m-%d %H:%M")

#calculate distance in hours of beginning and end of foraging from peak low tide
shells$tide_start <- difftime ( shells$start_time, shells$tide_time, units = "hours" )
shells$tide_end   <- difftime ( shells$end_time, shells$tide_time, units = "hours" )

#manually assign high tide value (based on averaged observed values for high tides)
high_tide <- 3 #in meters. 

#calculate average height of tide in trip with matrix approximation 
for ( i in 1:nrow(shells)) {
  if ( !is.na(shells$tide_start[i])){
    seq_time <- seq(as.numeric(shells$tide_start[i]), as.numeric(shells$tide_end[i]), 0.1)
    depth_at_seq_time <- exp(-0.065 * (seq_time) ^2 ) * (shells$tide_height[i] - high_tide) + high_tide
    shells$avg_tide_depth[i] <- mean(depth_at_seq_time)
  } else {
    shells$avg_tide_depth[i] <- NA
  }
}

##########################
#Cleaning the What column#
##########################
remove <- data.frame( to_replace = c("fukuilile", "fukulike", 
                                     "kichamvi", "kichonvi", "vichomvi",
                                     "kijnu", "vijinu",
                                     "kome 12 korongonjo", "similkorongonjo",
                                     "small ngisi", "uma"),
                      replacement = c("fukulile", "fukulile",
                                      "kichomvi", "kichomvi", "kichomvi",
                                      "kijinu", "kijinu",
                                      "kome, 12 korongonjo", "korongonjo",
                                      "ngisi", "una"))

for ( i in 1:nrow(remove)) shells$what <- str_replace(shells$what, 
                                                      remove$to_replace[i], 
                                                      remove$replacement[i])

shell_items <- gsub('[[:digit:]]+', '', unique(shells$what))
shell_items <- strsplit (shell_items, ", ")
shell_items <- Reduce(c,shell_items)
shell_items <- unique(gsub(" ","",shell_items))
shell_items <- shell_items[-which(is.na(shell_items))]

for ( i in 1:nrow(shells)) {
  shells$n_item_types[i] <- sum (str_detect(shells$what[i], shell_items ))
}
shells$n_item_types[which(is.na(shells$n_item_types))] <- 0

#shells generated quantities

generated_quantities$n_shell_trips <- length(unique(shells$trip)) #number of shellfish foraging trips 
generated_quantities$n_shell_trips_person <- nrow(shells) #number of foraging trips/person 
generated_quantities$n_shell_ppl <- length(unique(shells$who))

######################
#CREATE LIST OF TRAPS----------------------------------------------------------------------------------------
######################
#two different datasets are prepared for traps:
#1- all_trap_data, which reports data for each observation of a trap - again who built it, but then how long since last visit to the trap. They are used with bernoulli models. These data are incomplete, though, as on days when no observer was present to record the observation of traps, the data is missing and unrecoverable
#2- traps, which reports data at the trap level -e.g. who built it, how long it was deployed etc- and is used in poisson models 
#note that during cleaning of the data, several problem arose, wich are carefully described in the non-shared corrections file.

####
#ISSUES ON TRAP DATA
####

#TRAPS BUILT AND NEVER CHECKED
#25 traps are recorded only on the day of construction. Often it is a last trap of a series, or in one case a whole full series. These are likely been forgotten and not checked again. They have zero returns

#CAPTURES FOR WHICH WE DO NOT HAVE INFO ON WHO BUILT THE TRAP
#the return data include three preys for which we don't know the trap 
#trapReturns[ which(is.na( trapReturns$trap_ID )), ]
#Two traps that captured preys are not associated to people who built them (one is present on sheet, but no person is reported, the other is not reported, last of a series). Both are from trips when all traps were set by the same person 

#THE CONSTRUCTION OF SOME TRAPS HAS NOT BEEN RECORDED AND WAS ADDED AFTERWARDS
#trap 4BX4Y_7 does not appear in the list of traps built because there was no observer on that day, they were added in the correction stage
#traps 9PK9Y have been built during trip 9PK9M and their construction was not reported because there was no observer on that trip. Their construction was added in the correction stage
#traps built during trip 672WN were recorded during a following trip, including who mounted them
#traps built during trip 238TI were build by only one individual
#data for traps built during trips 303BJ, 2294B, 7CR4K, 6DF4R, 4MK2R, 4BX4Y, 6K9DA, 9C4MZ, 6M4LT, 6B9LB, 4K6KG, 9Z6KC, 6F9YZ, 7K9ZT were recorded even though no observer was present
#traps built during trip 4XK2Q were not recorded and do not appear afterwards
#Two trips were named 6DF4R on 17th of July and August. The first trip was not observed, but in the following trips, traps named 6DF4R were observed Construction of traps 6DF4S 1 to 4 was added a posteriori and these traps and relative trip was renamed 6DF4S
#d_trips[which(d_trips$observer == "none" & str_detect(d_trips$trip_type, "build" )),"trip"]


##########################
#1-WRANGLING TO CREATE DATA FOR INDIVIDUAL OBSERVATIONS OF TRAPS
##########################
#create a df with all activities concerning traps
trap_codes <- c ( "FY_BU", "FY_EM", "FY_CK", "FY_DI" )
all_trap_data <- d_activity %>% filter ( what %in% trap_codes )%>% select(entry_ID, trip, what, time, details)
all_trap_data$trap_ID <- str_extract ( all_trap_data$details, ".{5}_[:digit:]{1,}")#.{5} means the five characters preceding _, [:digit:]{1,} means one or more numbers after _

####
###CLEAN TRAP RECORD
####
#removes the entries which are marked by a trap code but do not report a trap ID
#note that a prey was captured in an unmarked trap (429QY_)
all_trap_data <- all_trap_data [ -which ( is.na ( all_trap_data$trap_ID )), ]
#remove three traps that are shared. 
all_trap_data <-  all_trap_data %>% filter ( !entry_ID %in% c("as_8","as_14","as_368") ) 
#remove entries for finishing time building traps
all_trap_data <-  all_trap_data %>% filter ( str_detect(all_trap_data$details, "finish|modify"  , negate = TRUE) )

#mark traps that captured non edible things as check
#all_trap_data[which (all_trap_data$what == "FY_EM"),]
all_trap_data[which (str_detect(all_trap_data$details, "komba|tumbili|remains")),"what"] <- "FY_CK"
#mark traps that captured edible things as empty trap
all_trap_data[which (str_detect(all_trap_data$details, "remove prey from trap|kwarara")),"what"] <- "FY_EM"
#correct remount data that are marked FY_BU
all_trap_data[which(all_trap_data$what == "FY_BU" & str_detect(all_trap_data$details, "remount")),"what"] <- "FY_CK"

#selected entries where a trap has been reported twice, once with an important entry 
to_remove <- c("as_678", "as_2275", "as_2497", "as_6170") 
all_trap_data <-  all_trap_data[ - which (all_trap_data$entry_ID %in% to_remove), ]

#removes entries where a trap has been marked twice in the same trip for the same thing (e.g. twice marked as checked)
all_trap_data <- all_trap_data %>% 
  distinct(trip, what, trap_ID, .keep_all = TRUE)

#Change activity for traps that were broken and then likely fixed
to_change <- c("as_2340", "as_4575", "as_6359", "as_6381")
all_trap_data$what[which (all_trap_data$entry_ID %in% to_change)] <- "FY_CK"
#traps 222TW_5, 666BT_8 429QY_7,303BJ_3, 7K9ZT_7 and 9KK9B_10 broke but were then likely fixed and used again


####
###ADD TIME AND DAY
####
#any trap trip not appearing in the trip dataframe?
#any( ! all_trap_data$trip %in% d_trips$trip) 
#add date
for (i in 1:nrow(all_trap_data)) {
  all_trap_data$date[i] <- d_trips$date[ which ( d_trips$trip == all_trap_data$trip[i])] 
}

#make column with date and time
all_trap_data$time_date <- all_trap_data$time

#input last written time for data points without time
#length(which(all_trap_data$time_date == ""))#3285
while(length(ind <- which(all_trap_data$time_date == "")) > 0){
  all_trap_data$time_date[ind] <- ifelse(all_trap_data$trip[ind] == all_trap_data$trip[ind-1],
                                         all_trap_data$time_date[ind -1],
                                         NA)
}

#input time of arrival at destination for data points still without time
#length(which(is.na(all_trap_data$time_date)))
na_time <- which(is.na(all_trap_data$time_date))
for (i in 1:length(na_time)){
  if(sum(d_activity$what == "AR_DE" & d_activity$trip == all_trap_data$trip[na_time[i]]) == 1){
    if(       d_activity$time[which(d_activity$what == "AR_DE" & d_activity$trip == all_trap_data$trip[na_time[i]])] != "" &
              !is.na(d_activity$time[which(d_activity$what == "AR_DE" & d_activity$trip == all_trap_data$trip[na_time[i]])] )){
      all_trap_data$time_date[na_time[i]] <- d_activity$time[
        which(d_activity$what == "AR_DE" & d_activity$trip == all_trap_data$trip[na_time[i]])]
    }#if
  }#if
}#i

#input inferred time for traps 
all_trap_data$time_date[which(all_trap_data$what == "FY_BU" & str_detect(all_trap_data$trap_ID, "491QS"))] <- "10:00" #traps built in the morning - no time
all_trap_data$time_date[which(all_trap_data$what == "FY_BU" & str_detect(all_trap_data$trap_ID, "157FE")&is.na(all_trap_data$time_date))] <- "10:00" #traps built in the morning - no time
all_trap_data$time_date[which( str_detect(all_trap_data$trip, "154JR"))] <- "9:18" #closest recorded time
all_trap_data$time_date[which( str_detect(all_trap_data$trip, "441TL")&is.na(all_trap_data$time_date))] <- "12:39" #closest recorded time
all_trap_data$time_date[which( str_detect(all_trap_data$trip, "7KC4M")&is.na(all_trap_data$time_date))] <- "15:33" #closest recorded time
all_trap_data$time_date[which( str_detect(all_trap_data$trip, "5M3YA")&is.na(all_trap_data$time_date))] <- "17:03" #closest recorded time
all_trap_data$time_date[which( str_detect(all_trap_data$trip, "3T4XM")&is.na(all_trap_data$time_date))] <- "16:00" #closest recorded time
all_trap_data$time_date[which( str_detect(all_trap_data$trip, "2C4WK")&is.na(all_trap_data$time_date))] <- "13:40" #closest recorded time

#enter time of the closest following entry from the same trip
while(length(ind <- which(is.na(all_trap_data$time_date))) > 0){
  all_trap_data$time_date[ind] <- ifelse(all_trap_data$trip[ind] == all_trap_data$trip[ind+1],
                                         all_trap_data$time_date[ind +1],
                                         "CHECK")
}
#no more entries without time
#any(all_trap_data$time_date == "CHECK")


#join time and date
all_trap_data$time_date <- ifelse(!is.na(all_trap_data$time_date),
                                  paste(all_trap_data$date, " ", all_trap_data$time_date, sep = ""),
                                  NA)



all_trap_data <- all_trap_data[order(all_trap_data$date),]

#add day number from beginning of data collection
all_trap_data$day <- difftime (all_trap_data$date, min(all_trap_data$date), units = "day") + 1

#calculate number of days since last check
all_trap_data$days_since_check <- NA
for (i in 1:length(unique(all_trap_data$trap_ID))) {
  trap_temp <- all_trap_data %>% filter(trap_ID == unique(all_trap_data$trap_ID)[i] )
  #trap_temp <- trap_temp %>% distinct(trip, .keep_all = TRUE)#removes duplicated values by same trip and trap
  for ( j in 1:nrow( trap_temp)){
    trap_temp$days_since_check[j] <- ifelse( trap_temp$what[j] == "FY_BU", NA, 
                                             trap_temp$day[j] - trap_temp$day[j-1])
    all_trap_data$days_since_check [ which(all_trap_data$entry_ID == trap_temp$entry_ID[j] )] <- trap_temp$days_since_check[j]
  }
}

#calculate hours since last check
all_trap_data$hours_since_check <- NA
for (i in 1:length(unique(all_trap_data$trap_ID))) {
  trap_temp <- all_trap_data %>% filter(trap_ID == unique(all_trap_data$trap_ID)[i] )
  for ( j in 1:nrow( trap_temp)){
    trap_temp$hours_since_check[j] <- ifelse( trap_temp$what[j] == "FY_BU", NA, 
                                              difftime( trap_temp$time_date[j], trap_temp$time_date[j-1], units = "hours"))
    all_trap_data$hours_since_check [ which(all_trap_data$entry_ID == trap_temp$entry_ID[j] )] <- trap_temp$hours_since_check[j]
  }
}

# Check negative and times below 1h
#all_trap_data[which(all_trap_data$days_since_check == 0),]
to_remove <- c("as_634", "as_632", "as_673", "as_874", "as_2224", "as_2271", "as_2353", "as_2883", "as_2905", "as_2953", "as_3532", "as_3835", "as_3979", "as_4012", "as_4285", "as_4286", "as_4523", "as_4541", "as_4647", "as_4668", "as_4852", "as_4958", "as_5651", "as_6201" )
all_trap_data <-  all_trap_data[ - which (all_trap_data$entry_ID %in% to_remove), ]
#trap 157FE_4 was stepped on and remounted, 
#the dog got captured in 157FE_12 which needed to be remounted, 
#on several days. two trips were undertaken generating observations at times closer than one day
#several traps were remounted after removing prey and remounts were removed from list


####
###ASSOCIATE TRAP CHECK WITH SUCCESS OF TRAP (by matching trip and trap and what)
####
all_trap_data$success <- NA
all_trap_data$returns <- NA

#make dataframe with returns from traps only
trapReturns <- d_returns[d_returns$trip_type == "fyuka_check",]
trapReturns <- trapReturns %>% filter(description %in% c("bird", "game"))

#add success
for(i in 1:nrow(all_trap_data)){
  if( all_trap_data[i, "what"]  == "FY_EM" &
      all_trap_data[i, "trip"] %in% trapReturns$trip &
      all_trap_data[i, "trap_ID"] %in% trapReturns$trap_ID) {
    all_trap_data$success[i] <- 1
  } else {
    all_trap_data$success[i] <- 0
  }
}

#add returns
for(i in 1:nrow(all_trap_data)){
  if( all_trap_data[i, "success"]  == 1) {
    all_trap_data$returns[i] <- trapReturns$weight_g [which(trapReturns$trip ==  all_trap_data[i, "trip"]&
                                                              trapReturns$trap_ID == all_trap_data[i, "trap_ID"])]
  } else {
    all_trap_data$returns[i] <- 0
  }
}

##########################
#2-WRANGLING TO CREATE DATA FOR TRAP LEVEL DATA
##########################

#Create list of built traps
traps <- d_activity %>% filter ( what == "FY_BU" )%>% select(entry_ID, trip, actor_ID, who, time, details)#add long lat if needed, but needs to be entered
#remove two traps that are shared. 
traps <-  traps %>% filter ( !entry_ID %in% c("as_8","as_14","as_368") ) 
#remove entries for finishing time building traps
traps <-  traps %>% filter ( str_detect(traps$details, "finish|remount|modify", negate = TRUE) )


#create column with trap IDs
traps$trap_ID <- str_extract(traps$details, ".{5}_[:digit:]{1,}")#.{5} means the five characters preceding _, [:digit:]{1,} means one or more numbers after _

#add success
for(i in 1:nrow(traps)){
  if( traps$trap_ID[i] %in% trapReturns$trap_ID) {
    traps$success[i] <- sum(trapReturns$trap_ID == traps$trap_ID[i], na.rm = TRUE)
  } else {
    traps$success[i] <- 0
  }
}

#add returns
for(i in 1:nrow(traps)){
  if( traps$trap_ID[i] %in% trapReturns$trap_ID) {
    traps$returns[i] <- sum( as.integer(trapReturns$weight_g [which (trapReturns$trap_ID == traps$trap_ID[i])]), na.rm = TRUE)
  } else {
    traps$returns[i] <- 0
  }
}

#add exposure - i.e. number of days the trap was deployed
for(i in 1:nrow(traps)){
  traps$build_time[i] <- min(all_trap_data$time_date[which(all_trap_data$trap_ID == traps$trap_ID[i])], na.rm = TRUE)
  traps$unmount_time[i] <- max(all_trap_data$time_date[which(all_trap_data$trap_ID == traps$trap_ID[i])], na.rm = TRUE)
  traps$exposure[i] <- difftime(traps$unmount_time[i], traps$build_time[i], units = "day")
}


#associate who built a trap to each row
all_trap_data$ID_who_built<-NA
for (i in 1:nrow(all_trap_data)) {
  if(all_trap_data$trap_ID[i] %in% traps$trap_ID){
    all_trap_data$ID_who_built[i] <- traps$actor_ID[ which(traps$trap_ID == all_trap_data$trap_ID[i]) ]
  }
}

####
###MAKE PLOT WHICH SHOWS LIVES OF TRAPS
####
#visual checks for lives of traps
png("../plots/traps.png", height = 16, width = 10, units = "cm", res = 500)
par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
plot( NULL, xlim = c( 1, max(all_trap_data$day) ),
      ylim = c(1, length(unique(all_trap_data$trap_ID))),
      xlab = "day", ylab = "trap index")

for ( i in 1: length(unique(all_trap_data$trap_ID))){
  trap <- all_trap_data %>% filter(trap_ID == unique(all_trap_data$trap_ID)[i] )
  lines( trap$day, rep( i, nrow(trap)), col = col.alpha("lightblue", 0.8))
  points( trap$day, rep( i, nrow(trap)), cex = 0.2,col = col.alpha("cornflowerblue", 0.2))
  if ( any (trap$what == "FY_BU" )) points(trap$day[trap$what == "FY_BU" ], i, pch = 16, cex = 0.2, col = col.alpha("#1482ac", 0.8))
  if ( any (trap$what == "FY_DI" )) points(trap$day[trap$what == "FY_DI" ], i, pch = 16, cex = 0.2, col = col.alpha("darkred", 0.8))
  if ( any (trap$what == "FY_EM" & trap$success == TRUE )) {
    points(trap$day[trap$what == "FY_EM" & trap$success == TRUE ],
           rep( i, length(trap$day[trap$what == "FY_EM" & trap$success == TRUE ])), pch = 16, cex = 0.3, col = col.alpha("orange", 0.8))}
  all_trap_data
}
dev.off()


##########
#generated quantities
##########

#n people who participated to trap trips
generated_quantities$n_trap_trip_participant <- length(d_participants$ID[which(d_participants$trip %in% all_trap_data$trip)])#total trips/participant
generated_quantities$n_trap_participants <- length(unique(d_participants$ID[which(d_participants$trip %in% all_trap_data$trip )]))#n unique participants

#n trips
trap_trips <- d_trips$trip[ which (str_detect(d_trips$trip_type, "fyuka"))]
generated_quantities$n_trap_trips <- length(trap_trips)


######################
#CREATE DATA FRAMES WITH INFO ABOUT PARTICIPANTS-------------------------------------------------------------
######################
#Shellfish
shell_ppl <- data.frame(ID = unique(shells$who),
                        age = rep(NA, length( unique(shells$who))),
                        sex = rep(NA, length( unique(shells$who))),
                        height = rep(NA, length( unique(shells$who))),
                        weight = rep(NA, length( unique(shells$who))),
                        grip = rep(NA, length( unique(shells$who))),
                        knowledge = rep(NA, length( unique(shells$who)),),
                        data = "shells"
)
for(i in 1:nrow(shell_ppl)){
  shell_ppl$age[i] <- census$age[ which (census$ID == shell_ppl$ID[i])]
}

for(i in 1:nrow(shell_ppl)){
  shell_ppl$sex[i] <- census$sex[ which (census$ID == shell_ppl$ID[i])]
}

#Traps
people_traps <- unique( na.omit( all_trap_data$ID_who_built ))
trap_ppl <- data.frame(ID = people_traps,
                       age = rep(NA, length( people_traps )),
                       sex = rep(NA, length( people_traps )),
                       height = rep(NA, length( people_traps )),
                       weight = rep(NA, length( people_traps )),
                       grip = rep(NA, length( people_traps )),
                       knowledge = rep(NA, length( people_traps )),
                       data = "traps"
)
for(i in 1:nrow(trap_ppl)){
  trap_ppl$age[i] <- census$age[ which (census$ID == trap_ppl$ID[i])]
}

for(i in 1:nrow(trap_ppl)){
  trap_ppl$sex[i] <- census$sex[ which (census$ID == trap_ppl$ID[i])]
}


#####
#generated quantities
#####
#proportion females
generated_quantities$prop_fem_shell <- sum(shell_ppl$sex == "f")/nrow(shell_ppl)
generated_quantities$prop_fem_trap  <- sum(trap_ppl$sex == "f")/nrow(trap_ppl)

#age range
generated_quantities$min_age_trap <- min(trap_ppl$age, na.rm = T)
generated_quantities$max_age_trap <- max(trap_ppl$age, na.rm = T)
generated_quantities$mean_age_trap <-mean(trap_ppl$age, na.rm = T)
generated_quantities$min_age_shell <- min(shell_ppl$age, na.rm = T)
generated_quantities$max_age_shell <- max(shell_ppl$age, na.rm = T)
generated_quantities$mean_age_shell <-mean(shell_ppl$age, na.rm = T)


generated_quantities$n_below19_shell <- sum(shell_ppl$age <=19, na.rm = T)
generated_quantities$n_above20_shell <- sum(shell_ppl$age >=20, na.rm = T)

#age ranges of people who participated to trap trips
participants_traps <- data.frame(ID = unique(d_participants$ID[which(d_participants$trip %in% trap_trips)]),
           age = rep(NA, length( unique(d_participants$ID[which(d_participants$trip %in% trap_trips)]) ))
)
for(i in 1:nrow(trap_ppl)){
  participants_traps$age[i] <- census$age[ which (census$ID == participants_traps$ID[i])]
}

generated_quantities$min_age_trap_participant <- min(participants_traps$age, na.rm = T)
generated_quantities$max_age_trap_participant <- max(participants_traps$age, na.rm = T)
generated_quantities$mean_age_trap_participant <-mean(participants_traps$age, na.rm = T)


#group sizes
#function to get mode of distribution
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

n_members_shells <- d_participants %>% filter(trip %in% shells$trip) %>% group_by(trip) %>% count()
generated_quantities$mean_n_participants_shells <- mean(n_members_shells$n)
generated_quantities$median_n_participants_shells <- median(n_members_shells$n)
generated_quantities$mode_n_participants_shells <- getmode(n_members_shells$n)

n_members_traps <- d_participants %>% filter(trip %in% trap_trips) %>% group_by(trip) %>% count()
generated_quantities$mean_n_participants_traps <- mean(n_members_traps$n)
generated_quantities$median_n_participants_traps <- median(n_members_traps$n)
generated_quantities$mode_n_participants_traps <- getmode(n_members_traps$n)




#####
#ADD ANTHROPOMETRICS
#####
##########################
#N anthropometric measures
##########################
d_anthropometrics$anonymeID <- anonyme(d_anthropometrics$ID)

count_nobs_anthro_shells <- d_anthropometrics %>% filter(anonymeID %in% shell_ppl$anonymeID) %>% count(anonymeID)
count_nobs_anthro_traps <- d_anthropometrics %>% filter(anonymeID %in% trap_ppl$anonymeID) %>% count(anonymeID)

generated_quantities$max_n_antrho_obs <- max( c( count_nobs_anthro_traps$n, count_nobs_anthro_shells$n))

#individuals were measured one to four times 

#add last measure per individual. 
anthropometrics <- d_anthropometrics %>% 
  mutate(date_measurement=as.Date(date_measurement, format= "%d/%m/%Y"))%>% 
  group_by(ID) %>%  
  arrange(desc(date_measurement)) %>%
  slice(1)
for(i in 1:nrow(anthropometrics)){
  anthropometrics$sex[i] <- census$sex[ which (census$ID == anthropometrics$ID[i])]
}

#####
#generated quantities
#####
#last measures were conducted within a time interval:
generated_quantities$latest_date_anthro <- max(anthropometrics$date_measurement)
generated_quantities$eaerliest_date_anthro <- min(anthropometrics$date_measurement)


###############
#ANONIMIZE ID
###############

shells$anonymeID <- anonyme(shells$who )
shell_ppl$anonymeID <- anonyme(shell_ppl$ID )
all_trap_data$anonymeID <- anonyme(all_trap_data$ID_who_built )
all_trap_data$anonymeID[ which(all_trap_data$anonymeID == "NA")] <- NA
traps$anonymeID <- anonyme(traps$actor_ID)
traps$anonymeID[ which(traps$anonymeID == "NA")] <- NA
trap_ppl$anonymeID <- anonyme(trap_ppl$ID )
anthropometrics$anonymeID <- anonyme(anthropometrics$ID)

###############
#make big dfs for all individuals with info
###############
#create df with all people from knowledge interview and anthropometric stuff, then merge to other dfs
knowledge_ppl <- data.frame( ID = rep(NA, d_knowledge$N),
                             age = d_knowledge$A,
                             sex = d_knowledge$S,
                             height = rep(NA, d_knowledge$N),
                             weight = rep(NA, d_knowledge$N),
                             grip = rep(NA, d_knowledge$N),
                             knowledge = rep(NA, d_knowledge$N),
                             data = rep("knowledge", d_knowledge$N),
                             anonymeID = rownames(d_knowledge$Y_l))

anthro_ppl <- data.frame( ID = anthropometrics$ID,
                          age = anthropometrics$age,
                          sex = anthropometrics$sex,
                          height = anthropometrics$height,
                          weight = anthropometrics$weight,
                          grip = anthropometrics$grip,
                          knowledge = rep(NA, nrow(anthropometrics)),
                          data = rep("anthropometrics", nrow(anthropometrics)),
                          anonymeID = anthropometrics$anonymeID)
anthro_ppl <- anthro_ppl[- which(is.na(anthro_ppl$age)),]

#bind and clean
shell_ppl <- rbind(shell_ppl, knowledge_ppl, anthro_ppl)
shell_ppl <- shell_ppl[!duplicated(shell_ppl$anonymeID), ]
trap_ppl <- rbind(trap_ppl, knowledge_ppl, anthro_ppl)
trap_ppl <- trap_ppl[!duplicated(trap_ppl$anonymeID), ]

#Shellfish
for(i in 1:nrow(shell_ppl)){
  if (shell_ppl$anonymeID[i] %in% anthropometrics$anonymeID){
    shell_ppl$height[i] <- anthropometrics$height[which(anthropometrics$anonymeID == shell_ppl$anonymeID[i])]
    shell_ppl$weight[i] <- anthropometrics$weight[which(anthropometrics$anonymeID == shell_ppl$anonymeID[i])]
    shell_ppl$grip[i] <- anthropometrics$grip[which(anthropometrics$anonymeID == shell_ppl$anonymeID[i])]
  } 
} 

#Traps
for(i in 1:nrow(trap_ppl)){
  if (trap_ppl$anonymeID[i] %in% anthropometrics$anonymeID){
    trap_ppl$height[i] <- anthropometrics$height[which(anthropometrics$anonymeID == trap_ppl$anonymeID[i])]
    trap_ppl$weight[i] <- anthropometrics$weight[which(anthropometrics$anonymeID == trap_ppl$anonymeID[i])]
    trap_ppl$grip[i] <- anthropometrics$grip[which(anthropometrics$anonymeID == trap_ppl$anonymeID[i])]
  } 
}



##################
#CONNECT KNOWLEDGE DATA
###############


list_knowledge <- data.frame(
  anonymeID = rownames(d_knowledge$Y_l),
  knowledge = rowSums(d_knowledge$Y_l)
)
for(i in 1:nrow(shell_ppl)){
  shell_ppl$knowledge[i] <- ifelse(shell_ppl$anonymeID[i] %in% list_knowledge$anonymeID, 
                                   list_knowledge$knowledge [which (list_knowledge$anonymeID == shell_ppl$anonymeID[i])], 
                                   NA)
}

for(i in 1:nrow(trap_ppl)){
  trap_ppl$knowledge[i] <- ifelse(trap_ppl$anonymeID[i] %in% list_knowledge$anonymeID, 
                                  list_knowledge$knowledge [which (list_knowledge$anonymeID == trap_ppl$anonymeID[i])], 
                                  NA)
}


#create new df with all answers for all ppl in the sample
shell_knowledge <- list(matrix(NA, nrow = nrow(shell_ppl), ncol = ncol(d_knowledge$Y_l), 
                          dimnames = list( shell_ppl$anonymeID, colnames(d_knowledge$Y_l)) ),
                        matrix(NA, nrow = nrow(shell_ppl), ncol = ncol(d_knowledge$Y_q), 
                               dimnames = list( shell_ppl$anonymeID, colnames(d_knowledge$Y_q)) ),
                        matrix(NA, nrow = nrow(shell_ppl), ncol = ncol(d_knowledge$Y_r), 
                               dimnames = list( shell_ppl$anonymeID, colnames(d_knowledge$Y_r)) ))
for(i in 1:nrow(shell_ppl)){
  if(shell_ppl$anonymeID[i] %in% rownames(d_knowledge$Y_l)){ 
    shell_knowledge[[1]][i,] <- d_knowledge$Y_l [which (rownames(d_knowledge$Y_l) == rownames(shell_knowledge[[1]])[i]),]}
}
for(i in 1:nrow(shell_ppl)){
  if(shell_ppl$anonymeID[i] %in% rownames(d_knowledge$Y_q)){ 
    shell_knowledge[[2]][i,] <- d_knowledge$Y_q [which (rownames(d_knowledge$Y_q) == rownames(shell_knowledge[[2]])[i]),]}
}
for(i in 1:nrow(shell_ppl)){
  if(shell_ppl$anonymeID[i] %in% rownames(d_knowledge$Y_r)){ 
    shell_knowledge[[3]][i,] <- d_knowledge$Y_r [which (rownames(d_knowledge$Y_r) == rownames(shell_knowledge[[3]])[i]),]}
}


trap_knowledge <- list(matrix(NA, nrow = nrow(trap_ppl), ncol = ncol(d_knowledge$Y_l), 
                               dimnames = list( trap_ppl$anonymeID, colnames(d_knowledge$Y_l)) ),
                        matrix(NA, nrow = nrow(trap_ppl), ncol = ncol(d_knowledge$Y_q), 
                               dimnames = list( trap_ppl$anonymeID, colnames(d_knowledge$Y_q)) ),
                        matrix(NA, nrow = nrow(trap_ppl), ncol = ncol(d_knowledge$Y_r), 
                               dimnames = list( trap_ppl$anonymeID, colnames(d_knowledge$Y_r)) ))
for(i in 1:nrow(trap_ppl)){
  if(trap_ppl$anonymeID[i] %in% rownames(d_knowledge$Y_l)){ 
    trap_knowledge[[1]][i,] <- d_knowledge$Y_l [which (rownames(d_knowledge$Y_l) == rownames(trap_knowledge[[1]])[i]),]}
}
for(i in 1:nrow(trap_ppl)){
  if(trap_ppl$anonymeID[i] %in% rownames(d_knowledge$Y_q)){ 
    trap_knowledge[[2]][i,] <- d_knowledge$Y_q [which (rownames(d_knowledge$Y_q) == rownames(trap_knowledge[[2]])[i]),]}
}
for(i in 1:nrow(trap_ppl)){
  if(trap_ppl$anonymeID[i] %in% rownames(d_knowledge$Y_r)){ 
    trap_knowledge[[3]][i,] <- d_knowledge$Y_r [which (rownames(d_knowledge$Y_r) == rownames(trap_knowledge[[3]])[i]),]}
}


######################
#PREPARE DATA FOR SAVING-------------------------------------------------------------------------------------
######################


#prepare table with shell data
tide_data <- shells %>% select(anonymeID, weight_g, trip_length, tide_height, avg_tide_depth, tide_time, tide_start, tide_end)

#prepare data TEMP
shells$success <- 1
shells <- shells %>% select(anonymeID, success, weight_g, trip_length, tide_height, avg_tide_depth, n_item_types)
traps <- traps %>% select(anonymeID, success, returns, exposure, trap_ID )
all_trap_data  <- all_trap_data  %>% select(anonymeID, success, returns, days_since_check, hours_since_check, trap_ID)
colnames(shells) <- c( "anonymeID", "success", "returns", "lenght_min", "tide_height_m", "tide_avg_depth", "n_item_types")
colnames(all_trap_data)  <- c( "anonymeID", "success", "returns", "lenght_day", "lenght_hour", "trap_ID")
shell_ppl <- shell_ppl %>% select(-"ID")
trap_ppl  <- trap_ppl %>% select(-"ID")



#input data for people who's age is missing - from consultation with research assistant
shell_ppl$age[which (shell_ppl$anonymeID == "12588")] <- 7
shell_ppl$age[which (shell_ppl$anonymeID == "252472")] <- 30
shell_ppl$age[which (shell_ppl$anonymeID == "84993")] <- 15
trap_ppl$age[which (trap_ppl$anonymeID == "117357")] <- 45

#####
#generated quantities
#####
generated_quantities$prop_complete_data_shells <- nrow(shells[ complete.cases(shells),])/nrow(shells) #60%
generated_quantities$prop_complete_data_traps <- nrow(traps[ complete.cases(traps),])/nrow(traps) #98%
generated_quantities$prop_complete_data_all_traps <- nrow(all_trap_data[ complete.cases(all_trap_data),])/nrow(all_trap_data) #84%, most of which are building times, which need to be removed anyway
generated_quantities$prop_miss_ppl_ID_traps <- sum(is.na(traps$anonymeID))/nrow(traps) #2%
generated_quantities$n_traps <- nrow(traps)#n traps
generated_quantities$n_trap_builders <- length(unique(traps$anonymeID))#n builders

#n traps captured something
generated_quantities$n_traps_success1 <- length(unique(traps$trap_ID[which(traps$success >= 1)]))#traps that captured
#total n success
generated_quantities$n_traps_preys <- sum(traps$success)




######################
#GENERATE DATA LIST AND SAVE---------------------------------------------------------------------------------
######################

d <- list (
  shells = shells[ complete.cases(shells),], 
  shell_ppl = shell_ppl,
  shell_k = shell_knowledge,
  all_traps = all_trap_data[ complete.cases(all_trap_data),],
  traps = traps[complete.cases(traps),],
  trap_ppl = trap_ppl,
  trap_k = trap_knowledge
)

rm(list=setdiff(ls(), c("d", "shells", "all_trap_data", "traps", "trap_ppl", "shell_ppl", "tide_data", "generated_quantities")))

list.save(d, '2_Data_preparation/processed_data.RData')
write.csv(tide_data, '2_data_preparation/tide_data.csv', row.names = FALSE)
list.save(generated_quantities, "2_data_preparation/generated_quantities.RData")

