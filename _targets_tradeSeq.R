##### setup #####
library(future)
library(targets)
library(magrittr)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)
options(future.globals.maxSize = 24000 * 1024^2)

source("./R/functions_tradeSeq.R")

tar_option_set(packages = c("qs", 
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
                            "scLANE", 
                            "gamlss", 
                            "tradeSeq", 
                            "magrittr", 
                            "parallel", 
                            "S4Vectors", 
                            "doParallel", 
                            "BiocParallel", 
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
# targets::tar_watch(level_separation = 1200, seconds = 120, seconds_max = 360, project = "tradeSeq")
# targets::tar_manifest(script = "_targets_tradeSeq.R") %>% 
# left_join((targets::tar_progress(store = "store_tradeSeq")), by = "name") %>% 
# mutate(progress = if_else(is.na(progress), "queued", progress)) %>% 
# count(progress) %>% 
# mutate(p = n / sum(n))

##### upstream targets #####
sims_single_subj <- data.frame(sim_file = list.files("store_simulation/objects/", pattern = "sim_single_*")) %>% 
                    dplyr::rowwise() %>% 
                    dplyr::mutate(ref_dataset = gsub("sim_single_", "", sim_file), 
                                  ref_dataset = gsub("_.*", "", ref_dataset), 
                                  dyn_gene_freq = gsub(paste0("sim_single_", ref_dataset, "_"), "", sim_file), 
                                  dyn_gene_freq = as.numeric(gsub("_.*", "", dyn_gene_freq)), 
                                  n_cells = as.numeric(gsub(paste0("sim_single_", ref_dataset, "_", dyn_gene_freq, "_"), "", sim_file)), 
                                  tradeseq_res_name = paste0("tradeSeq_", ref_dataset, "_DEG_", dyn_gene_freq, "_N_", n_cells)) %>% 
                    dplyr::ungroup()
sims_single_subj_symbol <- rlang::syms(sims_single_subj$sim_file)
sims_single_subj_file_symbol <- rlang::syms(paste0("file_", sims_single_subj$sim_file))
tradeSeq_symbols <- rlang::syms(sims_single_subj$tradeseq_res_name)

##### targets #####
list(
  tar_eval(values = list(symbol = sims_single_subj_file_symbol, 
                         file_string = paste0("store_simulation/objects/", sims_single_subj$sim_file)), 
           tar_target(symbol, 
                      file_string, 
                      format = "file", 
                      deployment = "main")), 
  tar_eval(values = list(symbol = sims_single_subj_symbol, 
                         file_symbol = sims_single_subj_file_symbol), 
           tar_target(symbol, qs::qread(file_symbol))), 
  tar_eval(values = list(symbol = tradeSeq_symbols, 
                         data_symbol = sims_single_subj_symbol), 
           tar_target(symbol, run_tradeSeq(data_symbol))), 
  tar_render(brain_metrics, "./Reports/tradeSeq_Brain_Metrics.Rmd"), 
  tar_render(endo_metrics, "./Reports/tradeSeq_Endocrinogenesis_Metrics.Rmd"), 
  tar_render(panc_metrics, "./Reports/tradeSeq_Pancreas_Metrics.Rmd")
)
