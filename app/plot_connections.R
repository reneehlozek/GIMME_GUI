# change file path to wherever testing folder is placed
setwd("XXX")

# reads in source code and stores function in your R workspace
source("location/to/your/gimmeGraph.R")

# only need to install once
#install.packages("R.matlab")
#install.packages("qgraph")
#install.packages("readxl")
#install.packages("igraph")
# check if installed
library(R.matlab)
library(qgraph)
library(readxl)
library(igraph)

# run function on input folder of xlsx matrices
# Need a "group" matrix even if empty
# Make sure FULL beta matrix is in input files - not just bottom half
gimmeGraph(input    = "betas",
           filetype = "csv",
           header   = TRUE,
           sep      = ",",
           weighted = TRUE)