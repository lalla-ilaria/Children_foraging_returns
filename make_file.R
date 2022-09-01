#load packages
library(rethinking)
library(dagitty)
library(rlist)
library(tidyverse)
library(ggnewscale)#for using two scales of color with ggplot
library(truncnorm)#for truncated normal
library(cowplot)#for pasting plots together
library(evoper)#for magnitude

#define colors and other plotting things
real_data <- list.load("2_data_preparation/processed_data.RData")
tide_data <- read.csv("2_data_preparation/tide_data.csv")
boycol  <- rgb(114/255, 181/255, 40/255) #"navyblue"
girlcol <- rgb(208/255, 29/255, 157/255) #"red3"
trapcol <- "#52b502" #"#5ca81e"
shellcol <-  "#ea5914"
othercol <- "grey30"##1482ac"
seq_trait <- seq(0,3,0.005)
age_plot <- 40


#compile models and create plots - simulated data
source("1_simulation/1_simulation.R")
source("1_simulation/2_analysis.R")


#compile models and create plots - real data
source("2_data_preparation/2_prep_data_functions.R")
source("3_analysis/1_analysis.R")
source("4_outcomes/1_plots.R")
source("4_outcomes/3_supp_info_plots.R")