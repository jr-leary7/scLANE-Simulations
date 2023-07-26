##### setup #####
library(future)
library(targets)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)
options(future.globals.maxSize = 24000 * 1024^2)

source("./R/functions_analysis.R")

tar_option_set(packages = c("qs", 
                            "stats", 
                            "igraph", 
                            "Matrix", 
                            "S4Vectors", 
                            "tidyverse", 
                            "rmarkdown", 
                            "tidymodels", 
                            "SummarizedExperiment", 
                            "SingleCellExperiment"), 
               error = "continue", 
               memory = "transient", 
               storage = "worker", 
               retrieval = "worker", 
               garbage_collection = TRUE, 
               seed = 312, 
               format = "qs")

##### monitoring #####
# targets::tar_watch(level_separation = 1200, seconds = 120, seconds_max = 360, project = "analysis")

##### targets #####
list(
  
)
