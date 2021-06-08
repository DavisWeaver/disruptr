#' function to compute null distribution of dnp
#'
#' `compute null` calculates a null distribution for the change in network potential for
#'  for each node in a cell signaling network.
#'
#'  The input for this function will be the output of [compute_dnp()].
#'   To compute the null distribution, the nodes in the provided cell signaling
#'   network will be randomly permuted `n` times, with dnp computed or each new
#'   cell signaling network. The mean and standard error of dnp for this set of random
#'   networks will constitute the null model that we will use for comparison.
#'   Be warned that this operation is extremely expensive computationally. It is
#'   recommended to either use a high-performance cluster or limit the computation of the
#'   null distribution to a small number of nodes.
#'   To distribute the workload over multiple cores, just specify ncores.
#'
#' @seealso [compute_dnp()] and [compute_np()]
#'
#' @param df output of [compute_dnp()]
#' @param n number of permutations
#' @param v character vector of protein names that we will compute the null
#'     distribution for, defaults to the top 50, as ranked by change in network potential from
#'     previous steps in the pipeline. Change n_genes parameter to compute for more or less
#' @param n_genes integer describing number of genes per sample that we will compute the null distribution for
#' @inheritParams compute_dnp
#'
#' @importFrom foreach %dopar% %:% %do%
#' @importFrom magrittr %>%
#'
#' @export

compute_null <- function(cache = NULL, df, ppi = "biogrid", n, all_genes = FALSE,
                         n_genes = 50, experiment_name, ncores = 1) {
  if(is.null(cache)) {
    stop("please provide a cache for saving and loading data/output")
  }
  #just return the file if we've already done this
  if(file.exists(paste0(cache, experiment_name, "dnpNull.Rda"))) {
    load(paste0(cache, experiment_name, "dnpNull.Rda"))
    return(df_np)
  }

  #load ppi
  g <- load_ppi(cache = cache, min_score = min_score, ppi = ppi)

  samples <- unique(df$sample_name)

  #grab top n genes per sample to iterate over
  v_list = get_topn(df =df, n_genes = n_genes)

  #get rid of dnp
  df <- df %>% dplyr::select(-dnp)
  #do we want to do a grouped split on sample and then iterate directly through the list? probably
  df_list <- df %>% dplyr::group_by(.data$sample_name) %>%
    dplyr::group_split()

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  out <-
    foreach::foreach(j = iterators::iter(df_list), .packages = "disruptr") %dopar% {
      cell_line_df <-
        foreach::foreach(i = 1:n, .combine = "rbind", .packages = "disruptr") %do%
        {
          #scramble the graph
          g_i <- get_random_graph(g)

          #calc change in network potential
          dnp_ij <- calc_dnp_i(j, g_i)
          dnp_ij$i = i
          return(dnp_ij)
        }

    }

  #close out parallel execution
  parallel::stopCluster(cl)



}

#' Helper function for compute_null - returns a graph with randomly permuted edges.
#'
#' currently just a wrapper for igraph::rewire but may add more functionality in the future
#'
#' @param g graph to be permuted
#'
#' @seealso [igraph::rewire()]

get_random_graph <- function(g) {
  g2 <- igraph::rewire(g, igraph::keeping_degseq(niter = igraph::vcount(g)*10))
  return(g2)
}

#' Helper function for compute_null - returns the top n genes by dnp for each sample
#'
#' @importFrom magrittr %>%

get_topn <- function(df, n_genes) {
  df %>% dplyr::group_by(.data$sample_name) %>%
    dplyr::slice_max(order_by = .data$dnp, n = n_genes) %>%
    dplyr::select(.data$gene_name, .data$sample_name) %>%
    dplyr::group_split()
}
