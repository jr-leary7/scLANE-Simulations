# run GLM scLANE on subsampled single-subject data
run_scLANE_GLM <- function(sim.data = NULL,
                           n.genes.sample = 2000,
                           n.cores = 4L) {
  # check inputs
  if (is.null(sim.data)) { stop("You failed to provide a SingleCellExperiment object to run_scLANE_GLM().") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }

  # prepare results list
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters",
                       "start_time",
                       "end_time",
                       "time_diff",
                       "mem_usage",
                       "testDynamic_results_raw",
                       "testDynamic_results_tidy",
                       "testSlope_results",
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
  mem_usage <- pryr::mem_change({
    start_time <- Sys.time()
    cell_offset <- createCellOffset(sim.data)
    gene_stats <- testDynamic(sim.data,
                              pt = pt_df,
                              genes = samp_genes,
                              size.factor.offset = cell_offset,
                              n.cores = n.cores,
                              approx.knot = TRUE,
                              verbose = FALSE)
    global_test_results <- getResultsDE(gene_stats,
                                        p.adj.method = "fdr",
                                        fdr.cutoff = 0.01) %>%
                           dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                              as.data.frame() %>%
                                              dplyr::mutate(gene = rownames(.))),
                                             by = c("Gene" = "gene"))
    end_time <- Sys.time()
    time_diff <- end_time - start_time
  })
  rmse_est <- purrr::imap(gene_stats, \(x, y) {
    rmse <- try({
      yardstick::rmse_vec(truth = as.numeric(BiocGenerics::counts(sim.data[y, ])), estimate = exp(x$Lineage_A$MARGE_Preds$marge_link_fit))
    }, silent = TRUE)
  })
  # set up results list & return
  res_list$sim_parameters <- list(n_cells = ncol(sim.data),
                                  n_genes = nrow(sim.data),
                                  n_dyn_genes = sum(SummarizedExperiment::rowData(sim.data)$geneStatus == "Dynamic"),
                                  p_dyn_genes = mean(SummarizedExperiment::rowData(sim.data)$geneStatus == "Dynamic"),
                                  method = "scLANE - GLM")
  res_list$start_time <- start_time
  res_list$end_time <- end_time
  res_list$time_diff <- time_diff
  res_list$mem_usage <- mem_usage
  res_list$testDynamic_results_raw <- gene_stats
  res_list$testDynamic_results_tidy <- global_test_results
  res_list$RMSE_estimates <- rmse_est
  return(res_list)
}
