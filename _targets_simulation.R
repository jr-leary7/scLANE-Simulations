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
                            "gtools", 
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
# targets::tar_manifest(script = "_targets_simulation.R") %>% 
# left_join((targets::tar_progress(store = "store_simulation")), by = "name") %>% 
# mutate(progress = if_else(is.na(progress), "queued", progress)) %>% 
# count(progress) %>% 
# mutate(p = n / sum(n))

##### simulation parameters #####

sim_params_single_subj <- tidyr::expand_grid(perc_dyn = c(1, 5, 10, 20), 
                                             n_cells = c(100, 250, 500, 1000, 3000, 5000))
sim_params_multi_subj <- tidyr::expand_grid(perc_dyn = c(5, 10), 
                                            n_cells = c(250, 500, 1000, 3000, 5000), 
                                            n_subjects = 6, 
                                            percent_overlap = c(70, 80, 90), 
                                            sample_alloc = c("balanced", "unbalanced"))
sim_params_multi_subj_het <- tidyr::expand_grid(perc_dyn = c(5, 10), 
                                                n_cells = c(250, 500, 1000, 3000, 5000), 
                                                n_subjects = 6, 
                                                percent_overlap = c(70, 80, 90), 
                                                sample_alloc = c("balanced", "unbalanced"), 
                                                group_overlap = c(40, 50, 60))
##### targets #####
list(
  # reference datasets 
  tar_target(brain_data_clean, fetch_lamanno_brain_data()), 
  tar_target(panc_data_clean, fetch_baron_pancreas_data()),
  tar_target(endo_data_clean, fetch_bastidas_ponce_pancreas_data()), 
  # single-subject simulations
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
                                                              n.cells = n_cells))), 
  # multi-subject simulations -- homogeneous trajectories 
  tar_map(values = sim_params_multi_subj, 
          tar_target(sim_multi_brain, simulate_multi_subject(ref.dataset = brain_data_clean, 
                                                             perc.dyn.genes = perc_dyn,
                                                             n.cells = n_cells, 
                                                             n.subjects = n_subjects, 
                                                             sample.alloc = sample_alloc, 
                                                             perc.overlap = percent_overlap))), 
  tar_map(values = sim_params_multi_subj, 
          tar_target(sim_multi_panc, simulate_multi_subject(ref.dataset = panc_data_clean, 
                                                            perc.dyn.genes = perc_dyn,
                                                            n.cells = n_cells, 
                                                            n.subjects = n_subjects, 
                                                            sample.alloc = sample_alloc, 
                                                            perc.overlap = percent_overlap))), 
  tar_map(values = sim_params_multi_subj, 
          tar_target(sim_multi_endo, simulate_multi_subject(ref.dataset = endo_data_clean, 
                                                            perc.dyn.genes = perc_dyn,
                                                            n.cells = n_cells, 
                                                            n.subjects = n_subjects, 
                                                            sample.alloc = sample_alloc, 
                                                            perc.overlap = percent_overlap))), 
  # multi-subject simulations -- heterogeneous trajectories 
  tar_map(values = sim_params_multi_subj_het, 
          tar_target(sim_multi_het_brain, simulate_multi_subject_het(ref.dataset = brain_data_clean, 
                                                                     perc.dyn.genes = perc_dyn,
                                                                     n.cells = n_cells, 
                                                                     n.subjects = n_subjects, 
                                                                     sample.alloc = sample_alloc, 
                                                                     perc.overlap = percent_overlap, 
                                                                     perc.overlap.group = group_overlap))), 
  tar_map(values = sim_params_multi_subj_het, 
          tar_target(sim_multi_het_panc, simulate_multi_subject_het(ref.dataset = panc_data_clean, 
                                                                    perc.dyn.genes = perc_dyn,
                                                                    n.cells = n_cells, 
                                                                    n.subjects = n_subjects, 
                                                                    sample.alloc = sample_alloc, 
                                                                    perc.overlap = percent_overlap, 
                                                                    perc.overlap.group = group_overlap))), 
  tar_map(values = sim_params_multi_subj_het, 
          tar_target(sim_multi_het_endo, simulate_multi_subject_het(ref.dataset = endo_data_clean, 
                                                                    perc.dyn.genes = perc_dyn,
                                                                    n.cells = n_cells, 
                                                                    n.subjects = n_subjects, 
                                                                    sample.alloc = sample_alloc, 
                                                                    perc.overlap = percent_overlap, 
                                                                    perc.overlap.group = group_overlap)))
  # QC report
  #tar_render(QC_report, "Reports/Simulation_QC.Rmd")
)
