get_pop_fit_Lamian <- function(lamian.res = NULL) {
  # check inputs
  if (is.null(lamian.res)) { stop("lamian.res must be non-NULL.") }

  # get predicted values from GAMs per-gene
  design <- lamian.res$design
  design <- design[, 1, drop = FALSE]
  knotnum <- lamian.res$knotnum
  pt_vec <- lamian.res$pseudotime
  pt_vec <- pt_vec[order(pt_vec)]
  gene <- rownames(lamian.res$statistics)
  fit_list <- purrr::map(gene, function(g) {
    tmp <- matrix(lamian.res$parameter[[g]]$beta, ncol = knotnum[g] + 4)
    beta <- as.vector(tmp[1, ])
    x <- purrr::map(rownames(design), function(i) {
      kronecker(diag(knotnum[g] + 4), design[i, , drop = FALSE])
    })
    if (knotnum[g] == 0) {
      phi <- splines::bs(pt_vec, intercept = TRUE)
    } else {
      knots <- seq(min(pt_vec), max(pt_vec), length.out = knotnum[g] + 2)[2:(knotnum[g] + 1)]
      phi <- splines::bs(pt_vec, knots = knots, intercept = TRUE)
    }
    i <- x[[1]]
    if (ncol(phi) == nrow(i)) {
      fit <- phi %*% i %*% beta
    } else {
      fit <- phi %*% t(i) %*% beta
    }
    fit_df <- data.frame(fit_link = as.numeric(fit))
    rownames(fit_df) <- rownames(fit)
    return(fit_df)
  })
  names(fit_list) <- gene
  return(fit_list)
}

get_preds_Lamian <- function(lamian.res = NULL) {
  # check inputs
  if (is.null(lamian.res)) { stop("lamian.res must be non-NULL.") }

  # run
  expr <- lamian.res$expr
  design <- lamian.res$design
  cellanno <- lamian.res$cellanno
  knotnum <- lamian.res$knotnum
  pt_vec <- lamian.res$pseudotime[colnames(expr)]
  gene <- rownames(expr)

  phi_list <- purrr::map(sort(unique(knotnum)), function(k) {
    if (k == 0) {
      phi <- splines::bs(pt_vec, intercept = TRUE)
    } else {
      knots <- seq(min(pt_vec), max(pt_vec), length.out = k + 2)[2:(k + 1)]
      phi <- splines::bs(pt_vec, knots = knots, intercept = TRUE)
    }
  })
  names(phi_list) <- as.character(sort(unique(knotnum)))

  sname <- purrr::map(rownames(design), function(i) {
    cellanno[cellanno[, 2] == i, 1]
  })
  names(sname) <- rownames(design)

  k <- unique(knotnum)[2]
  pred <- purrr::map(unique(knotnum), function(k) {
    genesub <- names(knotnum)[knotnum == k]
    B <- purrr::map(genesub, function(g) {
      lamian.res$parameter[[g]]$beta
    }) %>%
      purrr::reduce(rbind)
    if (length(genesub) == 1) {
      B <- matrix(B, nrow = 1)
    }
    rownames(B) <- genesub

    omega <- purrr::map(genesub, function(g) {
      lamian.res$parameter[[g]]$omega
    }) %>%
      purrr::reduce(rbind)
    if (length(genesub) == 1) {
      omega <- matrix(omega, nrow = 1)
    }
    rownames(omega) <- genesub

    phi <- phi_list[[as.character(k)]]
    phi <- purrr::map(rownames(design), function(ss) {
      phi[sname[[ss]], ]
    })
    names(phi) <- rownames(design)

    xs <- purrr::map(rownames(design), function(i) {
      kronecker(diag(k + 4), design[i, 1, drop = FALSE])
    })
    names(xs) <- rownames(design)

    phi_phi <- purrr::map(rownames(design), function(s) {
      t(phi[[s]]) %*% phi[[s]]
    })
    names(phi_phi) <- rownames(design)

    phiX <- purrr::map(rownames(design), function(s) {
      phi[[s]] %*% t(xs[[s]])
    })
    names(phiX) <- rownames(design)

    s <- rownames(design)[1]
    predtmp <- purrr::map(rownames(design), function(s) {
      sexpr <- expr[genesub, , drop = FALSE]
      sexpr_phibx <- sexpr[genesub, cellanno[, 2] == s, drop = FALSE] - (B[genesub, ] %*% t(phiX[[s]]))
      nb <- k + 4
      oinv <- purrr::map(genesub, function(g) {
        chol2inv(chol(matrix(omega[g, , drop = FALSE], nrow = nb)))
      })
      names(oinv) <- genesub
      Jchol <- purrr::map(genesub, function(g) {
        chol(phi_phi[[s]] + oinv[[g]])
      })
      names(Jchol) <- genesub
      Jsolve <- purrr::map(genesub, function(g) {
        as.numeric(chol2inv(Jchol[[g]]))
      }) %>%
        purrr::reduce(cbind)
      if (length(genesub) == 1) {
        Jsolve <- matrix(Jsolve, ncol = 1)
      }
      colnames(Jsolve) <- genesub

      K <- tcrossprod(t(phi[[s]]), sexpr_phibx)
      JK <- rowsum((Jsolve * K[rep(seq_len(nb), nb), , drop = FALSE]), rep(seq_len(nb), each =  nb))
      res <- t(phi[[s]] %*% JK)
      return(res)
    }) %>%
      purrr::reduce(cbind)
  }) %>%
    purrr::reduce(rbind)
  pred <- pred[gene, colnames(expr), drop = FALSE]
  populationFit <- t(get_pop_fit_Lamian(lamian.res = lamian.res) %>% purrr::reduce(cbind))
  pred_res <- pred + populationFit
  return(pred_res)
}

compute_resp_preds_Lamian <- function(lamian.preds, orig.sce) {
  # check inputs
  if (is.null(lamian.preds) || is.null(orig.sce)) { stop("Inputs must be non-NULL.") }
  size_factors <- BiocGenerics::sizeFactors(orig.sce)

  # compute fitted values on response scale (raw counts)
  pred_res_resp <- purrr::map(seq(nrow(lamian.preds)), function(r) {
    fitted_vals <- lamian.preds[r, , drop = FALSE]
    fitted_vals_resp <- t(t(2^(fitted_vals) - 1) * size_factors)
    return(fitted_vals_resp)
  }) %>%
    purrr::reduce(rbind)
  rownames(pred_res_resp) <- rownames(lamian.preds)
  return(pred_res_resp)
}

run_Lamian_GEE <- function(sim.data = NULL,
                           n.genes.sample = 2000,
                           n.cores = 4L,
                           lamian.iter = 100) {
  # check inputs
  if (is.null(sim.data)) { stop("You failed to provide necessary parameters to run_Lamian_GEE().") }
  if (lamian.iter <= 0) { stop("lamian.iter HAS to be positive, come on.") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }

  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters",
                       "start_time",
                       "end_time",
                       "time_diff",
                       "mem_usage",
                       "Lamian_results_raw",
                       "Lamian_results_tidy",
                       "RMSE_estimates")
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe -- {targets} takes care of seed
  p_dynamic <- SummarizedExperiment::rowData(sim.data) %>%
               as.data.frame() %>%
               dplyr::select(dplyr::contains("geneStatus_P")) %>%
               tidyr::pivot_longer(cols = tidyselect::everything(), values_to = "geneStatus") %>%
               dplyr::summarise(P = mean(geneStatus == "Dynamic")) %>%
               dplyr::pull(P)
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  if (!"geneDynamic_n" %in% colnames(SummarizedExperiment::rowData(sim.data))) {
    SummarizedExperiment::rowData(sim.data) <- SummarizedExperiment::rowData(sim.data) %>%
                                               as.data.frame() %>%
                                               dplyr::mutate(geneDynamic_n = rowSums(dplyr::across(tidyselect::contains("geneStatus_P"), \(x) x == "Dynamic"))) %>%
                                               S4Vectors::DataFrame()
  }
  samp_dyn_genes <- SummarizedExperiment::rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneDynamic_n > 0) %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- SummarizedExperiment::rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneDynamic_n == 0) %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]

  # prepare counts matrix & cell-ordering dataframe
  cell_anno <- SummarizedExperiment::colData(sim.data) %>%
               as.data.frame() %>%
               dplyr::mutate(Cell = rownames(.), .before = 1) %>%
               dplyr::select(Cell, Sample = subject_id) %>%
               dplyr::mutate(dplyr::across(tidyselect::everything(), as.character)) %>%
               magrittr::set_rownames(NULL)
  cell_pt <- as.integer(rank(sim.data$cell_time_normed))
  names(cell_pt) <- cell_anno$Cell
  samp_design <- data.frame(intercept = rep(1, length(unique(cell_anno$Sample)))) %>%
                 magrittr::set_rownames(unique(cell_anno$Sample)) %>%
                 as.matrix()
  sim_counts <- as.matrix(SingleCellExperiment::logcounts(sim.data))

  # run Lamian
    mem_usage <- pryr::mem_change({
      start_time <- Sys.time()
      gene_stats <- Lamian::lamian_test(expr = sim_counts,
                                        cellanno = cell_anno,
                                        pseudotime = cell_pt,
                                        design = samp_design,
                                        test.type = "time",
                                        test.method = "permutation",
                                        permuiter = lamian.iter,
                                        ncores = n.cores)
      global_test_results <- dplyr::rename(gene_stats$statistics,
                                           pvalue = pval.overall,
                                           Zstat = z.overall) %>%
                             dplyr::arrange(dplyr::desc(Zstat)) %>%
                             dplyr::mutate(pvalue_adj = stats::p.adjust(pvalue, method = "fdr")) %>%
                             dplyr::mutate(gene = rownames(.),
                                           gene_dynamic_overall = dplyr::if_else(pvalue_adj < 0.01, 1, 0)) %>%
                             dplyr::select(gene,
                                           Zstat,
                                           pvalue,
                                           pvalue_adj,
                                           gene_dynamic_overall) %>%
                             dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                                as.data.frame() %>%
                                                dplyr::select(contains("geneStatus_"), geneDynamic_n) %>%
                                                dplyr::mutate(gene = rownames(.))),
                                               by = "gene")
      end_time <- Sys.time()
      time_diff <- end_time - start_time
  })
  gene_stats$fitted_values_link <- get_preds_Lamian(lamian.res = gene_stats)
  gene_stats$fitted_values_resp <- compute_resp_preds_Lamian(lamian.preds = gene_stats$fitted_values_link, orig.sce = sim.data)
  rmse_est <- purrr::map(rownames(gene_stats$fitted_values_resp),
                         function(x) {
                           rmse <- try({
                             yardstick::rmse_vec(truth = sim_counts[x, ], estimate = gene_stats$fitted_values_resp[x, ])
                           }, silent = TRUE)
                         })
  names(rmse_est) <- rownames(gene_stats$fitted_values_resp)

  # set up results list & return
  res_list$sim_parameters <- list(n_cells = ncol(sim.data),
                                  n_genes = nrow(sim.data),
                                  n_dyn_genes = sum(SummarizedExperiment::rowData(sim.data)$geneDynamic_n > 0),
                                  p_dyn_genes = mean(SummarizedExperiment::rowData(sim.data)$geneDynamic_n > 0),
                                  method = "Lamian")
  res_list$start_time <- start_time
  res_list$end_time <- end_time
  res_list$time_diff <- time_diff
  res_list$mem_usage <- mem_usage
  res_list$Lamian_results_raw <- gene_stats
  res_list$Lamian_results_tidy <- global_test_results
  res_list$RMSE_estimates <- rmse_est
  return(res_list)
}
