##### setup #####
library(future)
library(targets)
library(magrittr)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)
options(future.globals.maxSize = 24000 * 1024^2)

source("./R/functions_scLANE_GLMM.R")

tar_option_set(packages = c("qs", 
                            "glm2", 
                            "mgcv", 
                            "pryr", 
                            "MASS", 
                            "rlang", 
                            "tidyr", 
                            "dplyr", 
                            "purrr", 
                            "stats", 
                            "broom", 
                            "scran", 
                            "gamlss", 
                            "scLANE", 
                            "glmmTMB", 
                            "foreach",  
                            "magrittr", 
                            "parallel", 
                            "bigstatsr", 
                            "S4Vectors", 
                            "doParallel", 
                            "SummarizedExperiment", 
                            "SingleCellExperiment"), 
               error = "continue", 
               memory = "transient",
               retrieval = "worker", 
               storage = "worker", 
               deployment = "worker", 
               garbage_collection = TRUE, 
               format = "qs")

##### monitoring #####
# targets::tar_watch(level_separation = 1200, seconds = 120, seconds_max = 360, project = "scLANE_GLMM")
# targets::tar_manifest(script = "_targets_scLANE_GLMM.R") %>% 
# left_join((targets::tar_progress(store = "store_scLANE_GLMM")), by = "name") %>% 
# mutate(progress = if_else(is.na(progress), "queued", progress)) %>% 
# count(progress) %>% 
# mutate(p = n / sum(n))

##### upstream targets #####
sims_multi_subj <- data.frame(sim_file = list.files("store_simulation/objects/", pattern = "sim_multi_*")) %>% 
                   dplyr::rowwise() %>% 
                   dplyr::mutate(ref_dataset = gsub("sim_multi_het_", "", sim_file), 
                                 ref_dataset = gsub("sim_multi_", "", ref_dataset), 
                                 ref_dataset = gsub("_.*", "", ref_dataset), 
                                 dyn_gene_freq = gsub(paste0("sim_multi_het_", ref_dataset, "_"), "", sim_file), 
                                 dyn_gene_freq = gsub(paste0("sim_multi_", ref_dataset, "_"), "", dyn_gene_freq), 
                                 dyn_gene_freq = as.numeric(gsub("_.*", "", dyn_gene_freq)), 
                                 n_cells = gsub(paste0("sim_multi_het_", ref_dataset, "_", dyn_gene_freq, "_"), "", sim_file), 
                                 n_cells = gsub(paste0("sim_multi_", ref_dataset, "_", dyn_gene_freq, "_"), "", n_cells), 
                                 n_cells = as.numeric(gsub("_.*", "", n_cells)), 
                                 n_subjects = gsub(paste0("sim_multi_het_", ref_dataset, "_", dyn_gene_freq, "_", n_cells, "_"), "", sim_file), 
                                 n_subjects = gsub(paste0("sim_multi_", ref_dataset, "_", dyn_gene_freq, "_", n_cells, "_"), "", n_subjects), 
                                 n_subjects = as.numeric(gsub("_.*", "", n_subjects)), 
                                 perc_overlap = gsub(paste0("sim_multi_het_", ref_dataset, "_", dyn_gene_freq, "_", n_cells, "_", n_subjects, "_"), "", sim_file), 
                                 perc_overlap = gsub(paste0("sim_multi_", ref_dataset, "_", dyn_gene_freq, "_", n_cells, "_", n_subjects, "_"), "", perc_overlap), 
                                 perc_overlap = as.numeric(gsub("_.*", "", perc_overlap)), 
                                 sample_alloc = gsub(paste0("sim_multi_het_", ref_dataset, "_", dyn_gene_freq, "_", n_cells, "_", n_subjects, "_", perc_overlap, "_"), "", sim_file), 
                                 sample_alloc = gsub(paste0("sim_multi_", ref_dataset, "_", dyn_gene_freq, "_", n_cells, "_", n_subjects, "_", perc_overlap, "_"), "", sample_alloc), 
                                 sample_alloc = gsub("_.*", "", sample_alloc), 
                                 group_overlap = dplyr::case_when(!grepl("_het_", sim_file) ~ NA_real_, 
                                                                  TRUE ~ as.numeric(gsub(".*_", "", sim_file))), 
                                 sclane_res_name = paste0("scLANE_GLMM_", ref_dataset, "_DEG_", dyn_gene_freq, "_N_", n_cells, "_SUBJ_", n_subjects, "_OVERLAP_", perc_overlap, "_", sample_alloc), 
                                 sclane_res_name = dplyr::case_when(!grepl("_het_", sim_file) ~ sclane_res_name, 
                                                                    TRUE ~ paste0(sclane_res_name, "_GROUP_", group_overlap))) %>% 
                   dplyr::ungroup() %>% 
                   dplyr::filter(n_cells <= 1000)
sims_multi_subj_symbol <- rlang::syms(sims_multi_subj$sim_file)
sims_multi_subj_file_symbol <- rlang::syms(paste0("file_", sims_multi_subj$sim_file))
scLANE_GLMM_symbols <- rlang::syms(sims_multi_subj$sclane_res_name)

##### targets #####
list(
  tar_eval(values = list(symbol = sims_multi_subj_file_symbol, 
                         file_string = paste0("store_simulation/objects/", sims_multi_subj$sim_file)), 
           tar_target(symbol, 
                      file_string, 
                      format = "file", 
                      deployment = "main")), 
  tar_eval(values = list(symbol = sims_multi_subj_symbol, 
                         file_symbol = sims_multi_subj_file_symbol), 
           tar_target(symbol, qs::qread(file_symbol))), 
  tar_eval(values = list(symbol = scLANE_GLMM_symbols, 
                         data_symbol = sims_multi_subj_symbol), 
           tar_target(symbol, run_scLANE_GLMM(data_symbol)))
)
