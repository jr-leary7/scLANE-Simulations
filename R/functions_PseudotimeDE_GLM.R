# run PseudotimeDE on subsampled single-subject data
run_PseudotimeDE <- function(sim.data = NULL,
                             n.genes.sample = 2000,
                             n.cores = 24L) {
  # check inputs
  if (is.null(sim.data)) { stop("You failed to provide necessary parameters to run_PseudotimeDE().") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  # prepare results list
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters",
                       "start_time",
                       "end_time",
                       "time_diff",
                       "mem_usage",
                       "PseudotimeDE_results_raw",
                       "PseudotimeDE_results_tidy",
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
  # run PsedotimeDE
  mem_usage <- pryr::mem_change({
    start_time <- Sys.time()
    bioc_par <- BiocParallel::MulticoreParam(workers = n.cores, RNGseed = 312)
    # run Slingshot to get a pseudotime for this dataset in PCA space
    rd <- irlba::prcomp_irlba(t(as.matrix(SingleCellExperiment::logcounts(sim.data))),
                              center = TRUE,
                              scale. = FALSE)
    sling_res <- slingshot::slingshot(rd$x[, 1:2])
    sling_pt <- slingshot::slingPseudotime(sling_res) %>%
                as.data.frame() %>%
                dplyr::mutate(cell = colnames(sim.data),
                              dplyr::across(dplyr::contains("Lineage"), \(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))),
                              .before = 1) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(pseudotime = mean(dplyr::c_across(dplyr::contains("Lineage")), na.rm = TRUE)) %>%
                dplyr::ungroup() %>%
                dplyr::select(-dplyr::contains("Lineage")) %>%
                dplyr::as_tibble()
    # generate 100 random subsamples of 80% of all cells
    subsample_cell_index <- BiocParallel::bplapply(seq(100), \(x) {
      set.seed(x)
      sampled_cells <- sample(x = seq(ncol(sim.data)),
                              size = 0.8 * ncol(sim.data),
                              replace = FALSE)
      sampled_cells
    }, BPPARAM = bioc_par)
    # re-run Slingshot on all subsamples
    subsample_sling_res <- BiocParallel::bplapply(seq(subsample_cell_index), \(x) {
      sce_sub <- sim.data[, subsample_cell_index[[x]]]
      rd_sub <- irlba::prcomp_irlba(t(as.matrix(SingleCellExperiment::logcounts(sce_sub))),
                                    center = TRUE,
                                    scale. = FALSE)
      sling_res_sub <- slingshot::slingshot(rd_sub$x[, 1:2])
      sling_pt_sub <- slingshot::slingPseudotime(sling_res_sub) %>%
                      as.data.frame() %>%
                      dplyr::mutate(cell = colnames(sce_sub),
                                    dplyr::across(dplyr::contains("Lineage"), \(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))),
                                    .before = 1) %>%
                      dplyr::rowwise() %>%
                      dplyr::mutate(pseudotime = mean(dplyr::c_across(dplyr::contains("Lineage")), na.rm = TRUE)) %>%
                      dplyr::ungroup() %>%
                      dplyr::select(-dplyr::contains("Lineage")) %>%
                      dplyr::as_tibble()
      merged_pt <- dplyr::left_join(sling_pt_sub,
                                    sling_pt,
                                    by = "cell")

      if (stats::cor(merged_pt$pseudotime.x, merged_pt$pseudotime.y) < 0) {
        sling_pt_sub <- dplyr::mutate(sling_pt_sub, pseudotime = 1 - pseudotime)
      }
      sling_pt_sub
    }, BPPARAM = bioc_par)
    BiocParallel::bpstop(bioc_par)
    # run PseudotimeDE
    gene_stats <- PseudotimeDE::runPseudotimeDE(samp_genes,
                                                ori.tbl = sling_pt,
                                                sub.tbl = subsample_sling_res,
                                                mat = as.matrix(BiocGenerics::counts(sim.data)),
                                                model = "nb",
                                                seed = 312,
                                                mc.cores = n.cores)
    # generate table of test statistics, adj. p-values, etc.
    global_test_results <-  dplyr::arrange(gene_stats, para.pv) %>%
                            dplyr::select(-gam.fit) %>%
                            dplyr::mutate(gene = samp_genes,
                                          pvalue_adj = stats::p.adjust(para.pv, method = "fdr"),
                                          gene_dynamic_overall = dplyr::if_else(pvalue_adj < 0.01, 1, 0)) %>%
                            dplyr::relocate(gene) %>%
                            dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                               as.data.frame() %>%
                                               dplyr::mutate(gene = rownames(.))),
                                              by = "gene")
    end_time <- Sys.time()
    time_diff <- end_time - start_time
  })
  rmse_est <- purrr::map(seq(nrow(gene_stats)), \(g) {
    rmse <- NA_real_  # maybe fix later but unimportant right now
  })
  names(rmse_est) <- samp_genes
  # set up results list & return
  res_list$sim_parameters <- list(n_cells = ncol(sim.data),
                                  n_genes = nrow(sim.data),
                                  n_dyn_genes = sum(SummarizedExperiment::rowData(sim.data)$geneStatus == "Dynamic"),
                                  p_dyn_genes = mean(SummarizedExperiment::rowData(sim.data)$geneStatus == "Dynamic"),
                                  method = "PseudotimeDE")
  res_list$start_time <- start_time
  res_list$end_time <- end_time
  res_list$time_diff <- time_diff
  res_list$mem_usage <- mem_usage
  res_list$PseudotimeDE_results_raw <- gene_stats
  res_list$PseudotimeDE_results_tidy <- global_test_results
  res_list$RMSE_estimates <- rmse_est
  return(res_list)
}
