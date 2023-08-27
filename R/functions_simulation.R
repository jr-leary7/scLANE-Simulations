# reference: https://www.rhondabacher.com/scaffold-vignette.pdf

# pull data from La Manno et al (2016) -- https://doi.org/10.1016/j.cell.2016.09.027
fetch_lamanno_brain_data <- function() {
  sce <- scRNAseq::LaMannoBrainData(which = "human-embryo")
  sce <- sce[rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, colSums(SingleCellExperiment::counts(sce)) > 0]
  sce <- sce[Rfast::Sort(rownames(sce)), ]  # necessary to play nice with scaffold
  return(sce)
}

# pull data from Baron et al (2016) -- https://doi.org/10.1016/j.cels.2016.08.011
fetch_baron_pancreas_data <- function() {
  sce <- scRNAseq::BaronPancreasData(which = "human")
  sce <- sce[rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, colSums(SingleCellExperiment::counts(sce)) > 0]
  sce <- sce[Rfast::Sort(rownames(sce)), ]  # necessary to play nice with scaffold
  return(sce)
}

# pull data from Bastidas-Ponce et al (2019) -- https://doi.org/10.1242/dev.173849
fetch_bastidas_ponce_pancreas_data <- function(conda_env_name = "/blue/rbacher/j.leary/py_envs/scLANE_sim_venv", conda_bin_path = "/apps/conda/23.3.1/condabin/conda") {
  # check inputs 
  if (!file.exists(conda_bin_path)) { stop("The provided conda binary doesn't exist.") }
  # fetch pancreatic endocrinogenesis data from scVelo
  reticulate::use_condaenv(condaenv = conda_env_name,
                           conda = conda_bin_path, 
                           required = TRUE)
  scvelo <- reticulate::import("scvelo")
  adata <- scvelo$datasets$pancreas()
  sce <- zellkonverter::AnnData2SCE(adata)
  # remove temporary download directory created by AnnData
  if (dir.exists("data/")) {
    system("rm -rf data/")
  }
  SingleCellExperiment::counts(sce) <- sce@assays@data$X
  sce@assays@data$X <- NULL
  sce@assays@data$spliced <- NULL
  sce@assays@data$unspliced <- NULL
  sce <- sce[Matrix::rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, Matrix::colSums(SingleCellExperiment::counts(sce)) > 0]
  sce <- sce[Rfast::Sort(rownames(sce)), ]  # necessary to play nice with scaffold
  return(sce)
}

# simulate data for one subject 
simulate_single_subject <- function(ref.dataset = NULL,
                                    perc.dyn.genes = NULL,
                                    n.cells = NULL, 
                                    n.knots = 2, 
                                    spline.degree = 2) {
  # check inputs
  if (is.null(ref.dataset) || is.null(perc.dyn.genes) || is.null(n.cells)) { stop("You're missing vital parameters for simulate_scaffold().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  perc.dyn.genes <- perc.dyn.genes * 0.01
  
  # set up simulation parameters
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  n_possible_dyn_genes <- ceiling((perc.dyn.genes / 0.8) * nrow(ref.dataset))  # make dynamic genes an 80% sample of the total pool of possible dynamic genes
  gene_means <- Matrix::rowMeans(SingleCellExperiment::counts(ref.dataset))
  Q80 <- unname(stats::quantile(gene_means, 0.8))
  high_exp_genes <- rownames(ref.dataset)[gene_means > Q80]
  if (length(high_exp_genes) < n_possible_dyn_genes) {
    Q70 <- unname(stats::quantile(gene_means, 0.7))
    high_exp_genes <- rownames(ref.dataset)[gene_means > Q70]
    if (length(high_exp_genes) < n_possible_dyn_genes) {
      stop("Your dataset has too few highly expressed genes to support the number of dynamic genes you want. Please reduce the % dynamic genes parameter.")
    }
  }
  # make sure pool of possible dynamic genes has high expression
  possible_dyn_genes <- sample(high_exp_genes,
                               size = n_possible_dyn_genes, 
                               replace = FALSE)
  dyn_genes <- sample(possible_dyn_genes, 
                      n_dyn_genes, 
                      replace = FALSE)
  # set up B-spline knots 
  my_knots <- purrr::map(seq(n_dyn_genes), \(k) {
    knot_1 <- stats::runif(1, 0.1, 0.5)
    knot_2 <- stats::runif(1, 0.5, 0.9)
    while (abs(knot_1 - knot_2) < 0.1) {
      knot_2 <- stats::runif(1, 0.5, 0.9)
    }
    other_knots <- c()
    if (n.knots > 2) {
      other_knots <- stats::runif(n.knots - 2, 0.1, 0.9)
    }
    knot_mat <- matrix(c(knot_1, knot_2, other_knots), 
                       nrow = 1, 
                       ncol = n.knots)
    return(knot_mat)
  })
  my_knots <- purrr::reduce(my_knots, rbind)
  rownames(my_knots) <- dyn_genes
  
  # set up B-spline theta (between-knot directions)
  ncol_theta <- n.knots + spline.degree + 1
  theta_pop <- stats::rnorm(5 * n_dyn_genes * ncol_theta, 5, 5)
  my_theta <- matrix(sample(theta_pop, ncol_theta * n_dyn_genes), 
                     ncol = ncol_theta, 
                     nrow = n_dyn_genes, 
                     byrow = TRUE)
  rownames(my_theta) <- dyn_genes
  dynamic_params <- list(propGenes = perc.dyn.genes,
                         dynGenes = dyn_genes, 
                         degree = spline.degree,
                         knots = my_knots,
                         theta = my_theta)
  
  # simulate 10X dataset
  scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.dataset,
                                                          sceUMI = TRUE,
                                                          useUMI = TRUE,
                                                          protocol = "droplet",
                                                          numCells = n.cells,
                                                          popHet = c(1, 1),
                                                          useDynamic = dynamic_params)
  sim_data <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
  SingleCellExperiment::counts(sim_data) <- Matrix::Matrix(sim_data@assays@data$umi_counts, sparse = TRUE)
  sim_data@assays@data$umi_counts <- NULL
  sim_data <- sim_data[rownames(ref.dataset), ]  # original ordering
  
  # typical scran + scater pre-processing pipeline
  sim_data <- sim_data[Matrix::rowSums(SingleCellExperiment::counts(sim_data) > 0) >= 3, Matrix::colSums(SingleCellExperiment::counts(sim_data)) > 0] 
  sim_data <- scater::logNormCounts(sim_data)
  var_decomp <- scran::modelGeneVar(sim_data)
  top2k_hvgs <- scran::getTopHVGs(var_decomp, n = 2000)
  sim_data <- scater::runPCA(sim_data,
                             subset_row = top2k_hvgs, 
                             ncomponents = 50)
  SingleCellExperiment::reducedDim(sim_data, "PCAsub") <- SingleCellExperiment::reducedDim(sim_data, "PCA")[, 1:30]
  sim_data <- scater::runUMAP(sim_data, 
                              dimred = "PCAsub", 
                              ncomponents = 2)
  g <- scran::buildSNNGraph(sim_data, 
                            use.dimred = "PCAsub", 
                            k = 30)
  clusters <- igraph::cluster_louvain(graph = g)$membership
  SingleCellExperiment::colLabels(sim_data) <- factor(clusters)
  SummarizedExperiment::colData(sim_data) <- SummarizedExperiment::colData(sim_data) %>%
                                             as.data.frame() %>%
                                             dplyr::mutate(cell_time = as.numeric(gsub("Cell_", "", rownames(.))),
                                                           cell_time_normed = (cell_time - min(cell_time)) / (max(cell_time) - min(cell_time)), 
                                                           subject = "P1") %>%
                                             S4Vectors::DataFrame()
  return(sim_data)
}

# simulate data for multiple subjects -- homogeneous trajectories
simulate_multi_subject <- function(ref.dataset = NULL,
                                   perc.dyn.genes = NULL,
                                   n.cells = NULL, 
                                   n.subjects = NULL, 
                                   sample.alloc = NULL, 
                                   perc.overlap = 75, 
                                   n.knots = 2, 
                                   spline.degree = 2) {
  # check inputs 
  if (is.null(ref.dataset) || is.null(perc.dyn.genes) || is.null(n.cells) || is.null(n.subjects) || is.null(sample.alloc)) { stop("You're missing vital parameters for simulate_multi_subject().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  
  # set up simulation parameters
  perc.dyn.genes <- perc.dyn.genes * 0.01
  perc.overlap <- perc.overlap * 0.01
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  n_possible_dyn_genes <- ceiling(n_dyn_genes * (1 / perc.overlap))
  
  # set up sample allocation -- balanced or unbalanced
  if (sample.alloc == "balanced") {
    n_per_subj <- rep(round(n.cells / n.subjects), n.subjects)
  } else {
    alloc_per_subj <- rep(0, n.subjects)
    while (any(alloc_per_subj < 0.1)) {
      alloc_per_subj <- as.numeric(gtools::rdirichlet(1, alpha = rep(10, n.subjects)))
    }
    n_per_subj <- round(alloc_per_subj * n.cells)
  }
  
  # make sure pool of possible dynamic genes has high expression
  gene_means <- Matrix::rowMeans(SingleCellExperiment::counts(ref.dataset))
  high_exp_quantile <- 0.8
  sample_Q <- unname(stats::quantile(gene_means, high_exp_quantile))
  high_exp_genes <- rownames(ref.dataset)[gene_means > sample_Q]
  while (length(high_exp_genes) < n_possible_dyn_genes && high_exp_quantile > 1e-5) {
    high_exp_quantile <- high_exp_quantile - 0.05
    sample_Q <- unname(stats::quantile(gene_means, high_exp_quantile))
    high_exp_genes <- rownames(ref.dataset)[gene_means > sample_Q]
  }
  if (length(high_exp_genes) < n_possible_dyn_genes) {
    stop("Your dataset has too few highly expressed genes to support the number of dynamic genes you want. Please reduce the % dynamic genes parameter.")
  }
  dyn_gene_pool <- sample(high_exp_genes, 
                          n_possible_dyn_genes, 
                          replace = FALSE)
  
  # set up B-spline knots 
  my_knots <- purrr::map(seq(n_possible_dyn_genes), \(k) {
    knot_1 <- stats::runif(1, 0.05, 0.5)
    knot_2 <- stats::runif(1, 0.5, 0.95)
    while (abs(knot_1 - knot_2) < 0.1) {
      knot_2 <- stats::runif(1, 0.5, 0.95)
    }
    knot_mat <- matrix(c(knot_1, knot_2), 
                       nrow = 1, 
                       ncol = 2)
    return(knot_mat)
  })
  my_knots <- purrr::reduce(my_knots, rbind)
  rownames(my_knots) <- dyn_gene_pool
  
  # set up B-spline theta (between-knot directions)
  ncol_theta <- n.knots + spline.degree + 1
  my_theta <- matrix(stats::rnorm(n_possible_dyn_genes * ncol_theta, 5, 5), 
                     ncol = ncol_theta, 
                     nrow = n_possible_dyn_genes, 
                     byrow = TRUE)
  rownames(my_theta) <- dyn_gene_pool
  
  # simulate 10X datasets per-subject
  subject_sims <- purrr::map(seq(n.subjects), \(s) {
    dyn_gene_sample <- sample(dyn_gene_pool, n_dyn_genes)
    dyn_gene_knots <- my_knots[dyn_gene_sample, ]
    dyn_gene_theta <- my_theta[dyn_gene_sample, ]
    dynamic_params <- list(propGenes = perc.dyn.genes,
                           dynGenes = dyn_gene_sample, 
                           degree = spline.degree,
                           knots = dyn_gene_knots,
                           theta = dyn_gene_theta)
    scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.dataset,
                                                            sceUMI = TRUE,
                                                            useUMI = TRUE,
                                                            protocol = "droplet",
                                                            numCells = n_per_subj[s],
                                                            popHet = c(1, 1),
                                                            useDynamic = dynamic_params)
    sim_data_subj <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
    SingleCellExperiment::counts(sim_data_subj) <- Matrix::Matrix(sim_data_subj@assays@data$umi_counts, sparse = TRUE)
    sim_data_subj@assays@data$umi_counts <- NULL
    sim_data_subj <- sim_data_subj[rownames(ref.dataset), ] 
    SummarizedExperiment::colData(sim_data_subj) <- SummarizedExperiment::colData(sim_data_subj) %>%
                                                    as.data.frame() %>%
                                                    dplyr::mutate(subject_id = paste0("P", s), 
                                                                  cell_time = as.numeric(gsub("Cell_", "", rownames(.))),
                                                                  cell_time_normed = (cell_time - min(cell_time)) / (max(cell_time) - min(cell_time))) %>%
                                                    S4Vectors::DataFrame()
    colnames(sim_data_subj) <- paste0("P", s, "_", colnames(sim_data_subj))
    colnames(SummarizedExperiment::rowData(sim_data_subj)) <- paste0("geneStatus_P", s)
    return(sim_data_subj)
  })
  
  # combine raw counts & metadata from multiple subjects 
  counts_mat <- purrr::map(subject_sims, counts) %>%
                purrr::reduce(Matrix::cbind2)
  cell_names <- purrr::map(subject_sims, colnames) %>% 
                purrr::reduce(c)
  gene_names <- rownames(subject_sims[[1]])
  rownames(counts_mat) <- gene_names
  colnames(counts_mat) <- cell_names
  row_data <- purrr::map(subject_sims, SummarizedExperiment::rowData) %>%
              purrr::reduce(cbind) %>%
              as.data.frame() %>% 
              dplyr::mutate(geneDynamic_n = rowSums(dplyr::across(tidyselect::contains("geneStatus_"), \(x) x == "Dynamic"))) %>% 
              S4Vectors::DataFrame()
  col_data <- purrr::map(subject_sims, SummarizedExperiment::colData) %>% 
              purrr::reduce(rbind)
  sim_data <- SingleCellExperiment::SingleCellExperiment(list(counts = counts_mat))
  SummarizedExperiment::colData(sim_data) <- col_data
  SummarizedExperiment::rowData(sim_data) <- row_data
  sim_data <- sim_data[rownames(ref.dataset), ] 
  
  # typical scran + scater pre-processing pipeline
  sim_data <- sim_data[Matrix::rowSums(SingleCellExperiment::counts(sim_data) > 0) >= 3, Matrix::colSums(SingleCellExperiment::counts(sim_data)) > 0] 
  sim_data <- scater::logNormCounts(sim_data)
  var_decomp <- scran::modelGeneVar(sim_data)
  top2k_hvgs <- scran::getTopHVGs(var_decomp, n = 2000)
  sim_data <- scater::runPCA(sim_data,
                             subset_row = top2k_hvgs, 
                             ncomponents = 50)
  SingleCellExperiment::reducedDim(sim_data, "PCAsub") <- SingleCellExperiment::reducedDim(sim_data, "PCA")[, 1:30]
  sim_data <- scater::runUMAP(sim_data, 
                              dimred = "PCAsub", 
                              ncomponents = 2)
  g <- scran::buildSNNGraph(sim_data, 
                            use.dimred = "PCAsub", 
                            k = 30)
  SingleCellExperiment::colLabels(sim_data) <- factor(igraph::cluster_louvain(g)$membership)
  return(sim_data)
}

# simulate data for multiple subjects -- heterogeneous trajectories
simulate_multi_subject_het <- function(ref.dataset = NULL,
                                       perc.dyn.genes = NULL,
                                       n.cells = NULL, 
                                       n.subjects = NULL, 
                                       sample.alloc = NULL, 
                                       perc.overlap = 75, 
                                       n.groups = 2,
                                       perc.overlap.group = 50, 
                                       n.knots = 2, 
                                       spline.degree = 2) {
  # check inputs 
  if (is.null(ref.dataset) || is.null(perc.dyn.genes) || is.null(n.cells) || is.null(n.subjects) || is.null(sample.alloc)) { stop("You're missing vital parameters for simulate_multi_subject().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  
  # set up simulation parameters
  subj_per_group <- n.subjects / n.groups
  perc.dyn.genes <- perc.dyn.genes * 0.01
  perc.overlap <- perc.overlap * 0.01
  perc.overlap.group <- perc.overlap.group * 0.01
  n_dyn_genes_per_subject <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  n_possible_dyn_genes_per_group <- ceiling(n_dyn_genes_per_subject * (1 / perc.overlap))
  n_possible_dyn_genes_total <- ceiling(n_possible_dyn_genes_per_group * (1 / perc.overlap.group))
  
  # set up sample allocation -- balanced or unbalanced
  if (sample.alloc == "balanced") {
    n_per_subj <- rep(round(n.cells / n.subjects), n.subjects)
  } else {
    alloc_per_subj <- rep(0, n.subjects)
    while (any(alloc_per_subj < 0.1)) {
      alloc_per_subj <- as.numeric(gtools::rdirichlet(1, alpha = rep(10, n.subjects)))
    }
    n_per_subj <- round(alloc_per_subj * n.cells)
  }
  n_per_subj <- split(n_per_subj, ceiling(seq_along(n_per_subj) / subj_per_group))

  # make sure pool of possible dynamic genes has high expression
  gene_means <- Matrix::rowMeans(SingleCellExperiment::counts(ref.dataset))
  high_exp_quantile <- 0.8
  sample_Q <- unname(stats::quantile(gene_means, high_exp_quantile))
  high_exp_genes <- rownames(ref.dataset)[gene_means > sample_Q]
  while (length(high_exp_genes) < n_possible_dyn_genes_total && high_exp_quantile > 1e-5) {
    high_exp_quantile <- high_exp_quantile - 0.05
    sample_Q <- unname(stats::quantile(gene_means, high_exp_quantile))
    high_exp_genes <- rownames(ref.dataset)[gene_means > sample_Q]
  }
  if (length(high_exp_genes) < n_possible_dyn_genes_total) {
    stop("Your dataset has too few highly expressed genes to support the number of dynamic genes you want. Please reduce the % dynamic genes parameter.")
  }
  dyn_gene_pool_total <- sample(high_exp_genes, 
                                n_possible_dyn_genes_total, 
                                replace = FALSE)
  
  # set up B-spline knots 
  my_knots <- purrr::map(seq(n_possible_dyn_genes_total), \(k) {
    knot_1 <- stats::runif(1, 0.05, 0.5)
    knot_2 <- stats::runif(1, 0.5, 0.95)
    while (abs(knot_1 - knot_2) < 0.1) {
      knot_2 <- stats::runif(1, 0.5, 0.95)
    }
    knot_mat <- matrix(c(knot_1, knot_2), 
                       nrow = 1, 
                       ncol = 2)
    return(knot_mat)
  })
  my_knots <- purrr::reduce(my_knots, rbind)
  rownames(my_knots) <- dyn_gene_pool_total
  
  # set up B-spline theta (between-knot directions)
  ncol_theta <- n.knots + spline.degree + 1
  my_theta <- matrix(stats::rnorm(n_possible_dyn_genes_total * ncol_theta, 5, 5), 
                     ncol = ncol_theta, 
                     nrow = n_possible_dyn_genes_total, 
                     byrow = TRUE)
  rownames(my_theta) <- dyn_gene_pool_total
  
  # create unique subject IDs
  subj_ids <- paste0("P", seq(n.subjects))
  subj_ids <- split(subj_ids, ceiling(seq_along(subj_ids) / subj_per_group))
  
  # simulate 10X dataset per-group, per-subject
  group_sims <- purrr::map(seq(n.groups), \(g) {
    dyn_gene_pool_group <- sample(dyn_gene_pool_total, n_possible_dyn_genes_per_group)
    subject_sims <- purrr::map(seq(n.subjects / n.groups), \(s) {
      subj_id <- subj_ids[[g]][s]
      dyn_gene_sample_subject <- sample(dyn_gene_pool_group, n_dyn_genes_per_subject)
      dyn_gene_subject_knots <- my_knots[dyn_gene_sample_subject, ]
      dyn_gene_subject_theta <- my_theta[dyn_gene_sample_subject, ]
      dynamic_params <- list(propGenes = perc.dyn.genes,
                             dynGenes = dyn_gene_sample_subject, 
                             degree = spline.degree,
                             knots = dyn_gene_subject_knots,
                             theta = dyn_gene_subject_theta)
      scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.dataset,
                                                              sceUMI = TRUE,
                                                              useUMI = TRUE,
                                                              protocol = "droplet",
                                                              numCells = n_per_subj[[g]][s],
                                                              popHet = c(1, 1),
                                                              useDynamic = dynamic_params)
      sim_data_subj <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
      SingleCellExperiment::counts(sim_data_subj) <- Matrix::Matrix(sim_data_subj@assays@data$umi_counts, sparse = TRUE)
      sim_data_subj@assays@data$umi_counts <- NULL
      sim_data_subj <- sim_data_subj[rownames(ref.dataset), ] 
      SummarizedExperiment::colData(sim_data_subj) <- SummarizedExperiment::colData(sim_data_subj) %>%
                                                      as.data.frame() %>%
                                                      dplyr::mutate(group_id = paste0("Group", g), 
                                                                    subject_id = subj_id, 
                                                                    cell_time = as.numeric(gsub("Cell_", "", rownames(.))),
                                                                    cell_time_normed = (cell_time - min(cell_time)) / (max(cell_time) - min(cell_time))) %>%
                                                      S4Vectors::DataFrame()
      colnames(sim_data_subj) <- paste0(subj_id, "_", colnames(sim_data_subj))
      colnames(SummarizedExperiment::rowData(sim_data_subj)) <- paste0("geneStatus_", subj_id)
      return(sim_data_subj)
    })
    return(subject_sims)
  })
  
  # combine raw counts & metadata from multiple subjects 
  counts_mat <- purrr::map(group_sims, \(g) {
    purrr::map(g, counts) %>% purrr::reduce(Matrix::cbind2)
  })
  counts_mat <- purrr::reduce(counts_mat, Matrix::cbind2)
  col_data <- purrr::map(group_sims, \(g) {
    purrr::map(g, SummarizedExperiment::colData) %>% purrr::reduce(rbind)
  })
  col_data <- purrr::reduce(col_data, rbind)
  row_data <- purrr::map(group_sims, \(g) {
    purrr::map(g, SummarizedExperiment::rowData) %>% 
      purrr::reduce(cbind) %>%
      as.data.frame() %>% 
      dplyr::mutate(!!paste0("geneDynamic_",  g[[1]]$group_id[1], "_n") := rowSums(dplyr::across(tidyselect::contains("geneStatus_"), \(x) x == "Dynamic")))
  })
  row_data <- purrr::reduce(row_data, cbind) %>% 
              dplyr::mutate(gene = rownames(.), .before = 1) %>% 
              S4Vectors::DataFrame()
  sim_data <- SingleCellExperiment::SingleCellExperiment(list(counts = counts_mat))
  SummarizedExperiment::colData(sim_data) <- col_data
  SummarizedExperiment::rowData(sim_data) <- row_data
  sim_data <- sim_data[rownames(ref.dataset), ] 
  
  # typical scran + scater pre-processing pipeline
  sim_data <- sim_data[Matrix::rowSums(SingleCellExperiment::counts(sim_data) > 0) >= 3, Matrix::colSums(SingleCellExperiment::counts(sim_data)) > 0] 
  sim_data <- scater::logNormCounts(sim_data)
  var_decomp <- scran::modelGeneVar(sim_data)
  top2k_hvgs <- scran::getTopHVGs(var_decomp, n = 2000)
  sim_data <- scater::runPCA(sim_data,
                             subset_row = top2k_hvgs, 
                             ncomponents = 50)
  SingleCellExperiment::reducedDim(sim_data, "PCAsub") <- SingleCellExperiment::reducedDim(sim_data, "PCA")[, 1:30]
  sim_data <- scater::runUMAP(sim_data, 
                              dimred = "PCAsub", 
                              ncomponents = 2)
  g <- scran::buildSNNGraph(sim_data, 
                            use.dimred = "PCAsub", 
                            k = 30)
  SingleCellExperiment::colLabels(sim_data) <- factor(igraph::cluster_louvain(g)$membership)
  return(sim_data)
}
