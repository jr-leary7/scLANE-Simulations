##### setup #####
library(targets)
library(tarchetypes)

source("./R/functions_scLANE_GLM.R")

tar_option_set(packages = c("qs", 
                            "glm2", 
                            "mgcv", 
                            "pryr", 
                            "MASS", 
                            "tidyr", 
                            "dplyr", 
                            "purrr", 
                            "stats", 
                            "broom", 
                            "scran", 
                            "gamlss", 
                            "scLANE", 
                            "foreach",  
                            "magrittr", 
                            "parallel", 
                            "bigstatsr", 
                            "S4Vectors", 
                            "doParallel", 
                            "SummarizedExperiment", 
                            "SingleCellExperiment"), 
               imports = c("scLANE"), 
               error = "continue", 
               memory = "transient",
               retrieval = "worker", 
               storage = "worker", 
               garbage_collection = TRUE, 
               format = "qs")

##### upstream targets #####
sims_single_subj <- daat.frame(sim_file = list.files("store_simulation/objects/", pattern = "sim_*"))

##### targets #####
list(
  tar_map(values = sims_single_subj, 
          tar_target(file, paste0("store_simulation/objects/", sim_file), format = "file")), 
  tar_map()
)
