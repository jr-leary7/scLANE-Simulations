# run tradeSeq on subsampled single-subject data
run_tradeSeq <- function(sim.data = NULL,
                         n.genes.sample = 2000,
                         n.cores = 4L) {
  # check inputs
  if (is.null(sim.data)) { stop("You failed to provide necessary parameters to run_tradeSeq().") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }

  # prepare results list
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters",
                       "start_time",
                       "end_time",
                       "time_diff",
                       "mem_usage",
                       "tradeSeq_results_raw",
                       "tradeSeq_results_tidy",
                       "RMSE_estimates")

  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe -- {targets} takes care of seed
  p_dynamic <- mean(SummarizedExperiment::rowData(sim.data)[, 1] == "Dynamic")
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- SummarizedExperiment::rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- SummarizedExperiment::rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  pt_df <- SummarizedExperiment::colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  sim_counts <- SingleCellExperiment::counts(sim.data[samp_genes, ])

  # run tradeSeq
  mem_usage <- pryr::mem_change({
    start_time <- Sys.time()
    bioc_par <- BiocParallel::MulticoreParam(workers = n.cores, RNGseed = 312)
    cell_offset <- scLANE::createCellOffset(sim.data)
    k_eval <- tradeSeq::evaluateK(sim_counts,
                                  pseudotime = pt_df,
                                  cellWeights = matrix(rep(1, nrow(pt_df)), ncol = 1),
                                  offset = log(1 / cell_offset),
                                  k = 3:10,
                                  plot = FALSE,
                                  nGenes = 500,
                                  verbose = FALSE,
                                  parallel = TRUE,
                                  BPPARAM = bioc_par)
    best_k <- c(3:10)[which.min(abs(colMeans(k_eval - rowMeans(k_eval))))]  # choose k w/ lowest MAD from mean AIC
    gene_stats <- tradeSeq::fitGAM(counts = sim_counts,
                                   pseudotime = pt_df,
                                   cellWeights = matrix(rep(1, nrow(pt_df)), ncol = 1),
                                   offset = log(1 / cell_offset),
                                   nknots = best_k,
                                   sce = FALSE,
                                   parallel = TRUE,
                                   verbose = FALSE,
                                   BPPARAM = bioc_par)
    BiocParallel::bpstop(bioc_par)
    global_test_results <- tradeSeq::associationTest(models = gene_stats, global = TRUE) %>%
                           dplyr::arrange(pvalue) %>%
                           dplyr::mutate(gene = rownames(.),
                                         pvalue_adj = stats::p.adjust(pvalue, method = "fdr"),
                                         gene_dynamic_overall = dplyr::case_when(pvalue_adj < 0.01 ~ 1, TRUE ~ 0)) %>%
                           dplyr::relocate(gene) %>%
                           dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                              as.data.frame() %>%
                                              dplyr::mutate(gene = rownames(.))),
                                             by = "gene")
    end_time <- Sys.time()
    time_diff <- end_time - start_time
  })
  rmse_est <- purrr::map(gene_stats, \(g) {
    rmse <- try({
      yardstick::rmse_vec(truth = g$y, estimate = g$fitted.values)
    }, silent = TRUE)
  })
  names(rmse_est) <- rownames(sim_counts)

  # set up results list & return
  res_list$sim_parameters <- list(n_cells = ncol(sim.data),
                                  n_genes = nrow(sim.data),
                                  n_dyn_genes = sum(SummarizedExperiment::rowData(sim.data)$geneStatus == "Dynamic"),
                                  p_dyn_genes = mean(SummarizedExperiment::rowData(sim.data)$geneStatus == "Dynamic"),
                                  method = "tradeSeq")
  res_list$start_time <- start_time
  res_list$end_time <- end_time
  res_list$time_diff <- time_diff
  res_list$mem_usage <- mem_usage
  res_list$tradeSeq_results_raw <- gene_stats
  res_list$tradeSeq_results_tidy <- global_test_results
  res_list$RMSE_estimates <- rmse_est
  return(res_list)
}
