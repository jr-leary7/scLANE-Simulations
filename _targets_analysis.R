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
                            "scLANE", 
                            "tradeSeq", 
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
  tar_target(file_Metrics_scLANE_GLM_Single_Subject_Endocrinogenesis, 
             "Data/Metrics_scLANE_GLM_Single_Subject_Endocrinogenesis", 
             format = "file", 
             deployment = "main"), 
  tar_target(file_Metrics_tradeSeq_Single_Subject_Endocrinogenesis, 
             "Data/Metrics_tradeSeq_Single_Subject_Endocrinogenesis", 
             format = "file", 
             deployment = "main"), 
  tar_target(file_Metrics_scLANE_GLM_Single_Subject_Pancreas, 
             "Data/Metrics_scLANE_GLM_Single_Subject_Pancreas", 
             format = "file", 
             deployment = "main"), 
  tar_target(file_Metrics_tradeSeq_Single_Subject_Pancreas, 
             "Data/Metrics_tradeSeq_Single_Subject_Pancreas", 
             format = "file", 
             deployment = "main"), 
  tar_target(file_Metrics_scLANE_GLM_Single_Subject_Brain, 
             "Data/Metrics_scLANE_GLM_Single_Subject_Brain", 
             format = "file", 
             deployment = "main"), 
  tar_target(file_Metrics_tradeSeq_Single_Subject_Brain, 
             "Data/Metrics_tradeSeq_Single_Subject_Brain", 
             format = "file", 
             deployment = "main"), 
  tar_target(Metrics_scLANE_GLM_Single_Subject_Endocrinogenesis, 
             qs::qread(file_Metrics_scLANE_GLM_Single_Subject_Endocrinogenesis)), 
  tar_target(Metrics_tradeSeq_Single_Subject_Endocrinogenesis, 
             qs::qread(file_Metrics_tradeSeq_Single_Subject_Endocrinogenesis)), 
  tar_target(Metrics_scLANE_GLM_Single_Subject_Pancreas, 
             qs::qread(file_Metrics_scLANE_GLM_Single_Subject_Pancreas)), 
  tar_target(Metrics_tradeSeq_Single_Subject_Pancreas, 
             qs::qread(file_Metrics_tradeSeq_Single_Subject_Pancreas)), 
  tar_target(Metrics_scLANE_GLM_Single_Subject_Brain, 
             qs::qread(file_Metrics_scLANE_GLM_Single_Subject_Brain)), 
  tar_target(Metrics_tradeSeq_Single_Subject_Brain, 
             qs::qread(file_Metrics_tradeSeq_Single_Subject_Brain)), 
  tar_target(metric_table_master, purrr::reduce(list(Metrics_scLANE_GLM_Single_Subject_Endocrinogenesis,
                                                     Metrics_tradeSeq_Single_Subject_Endocrinogenesis,
                                                     Metrics_scLANE_GLM_Single_Subject_Pancreas, 
                                                     Metrics_tradeSeq_Single_Subject_Pancreas, 
                                                     Metrics_scLANE_GLM_Single_Subject_Brain, 
                                                     Metrics_tradeSeq_Single_Subject_Brain), 
                                                rbind))
)
