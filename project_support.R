#load packages
library(rethinking)
library(dagitty)
library(rlist)
library(tidyverse)
library(ggnewscale)#for using two scales of color with ggplot
library(truncnorm)#for truncated normal
library(cowplot)#for pasting plots together
library(evoper)#for Magnitude()

library(tictoc)

#define colors and other plotting things
boycol  <- rgb(114/255, 181/255, 40/255) #"navyblue"
girlcol <- rgb(208/255, 29/255, 157/255) #"red3"
trapcol <- "#52b502" #"#5ca81e"
shellcol <-  "#ea5914"
othercol <- "grey30"##1482ac"
seq_trait <- seq(0,3,0.005)
age_plot <- 40

source("R/prep_data_functions.R")
source("R/simulation_functions.R")
