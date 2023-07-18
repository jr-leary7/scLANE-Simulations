##### setup #####
library(future)
library(targets)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)
options(future.globals.maxSize = 24000 * 1024^2)

source("./R/functions_simulation.R")

tar_option_set(packages = c("stats", 
                            "scran", 
                            "Rfast", 
                            "scater", 
                            "dplyr", 
                            "purrr", 
                            "igraph", 
                            "Matrix", 
                            "scRNAseq", 
                            "magrittr", 
                            "scaffold", 
                            "S4Vectors", 
                            "rmarkdown", 
                            "reticulate", 
                            "zellkonverter", 
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
# targets::tar_watch(level_separation = 1200, seconds = 120, seconds_max = 360, project = "simulation")

sim_params_single_subj <- tidyr::expand_grid(perc_dyn = c(1, 5, 10, 20), 
                                             n_cells = c(100, 250, 500, 1000, 5000))
##### targets #####
list(
  # reference datasets 
  tar_target(brain_data_clean, fetch_lamanno_brain_data()), 
  tar_target(panc_data_clean, fetch_baron_pancreas_data()),
  tar_target(endo_data_clean, fetch_bastidas_ponce_pancreas_data()), 
  # single subject simulations
  tar_map(values = sim_params_single_subj, 
          tar_target(sim_single_brain, simulate_single_subject(ref.dataset = brain_data_clean, 
                                                               perc.dyn.genes = perc_dyn,
                                                               n.cells = n_cells))), 
  tar_map(values = sim_params_single_subj, 
          tar_target(sim_single_panc, simulate_single_subject(ref.dataset = panc_data_clean, 
                                                              perc.dyn.genes = perc_dyn,
                                                              n.cells = n_cells))), 
  tar_map(values = sim_params_single_subj, 
          tar_target(sim_single_endo, simulate_single_subject(ref.dataset = endo_data_clean, 
                                                              perc.dyn.genes = perc_dyn,
                                                              n.cells = n_cells)))#, 
  # QC report
  #tar_render(QC_report, "Reports/Simulation_QC.Rmd")
)
