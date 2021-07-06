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
#' @param n_genes integer describing number of genes per sample that we will compute the null distribution for
#' @inheritParams compute_dnp
#'
#' @importFrom foreach %dopar% %do% %:%
#' @importFrom magrittr %>%
#'
#' @export

compute_null <- function(cache = NULL, df, ppi = "biogrid", n,
                         n_genes = 50, experiment_name, ncores = 4,
                         min_score = NULL) {
  if(is.null(cache)) {
    stop("please provide a cache for saving and loading data/output")
  }
  #just return the file if we've already done this
  if(file.exists(paste0(cache, experiment_name, "dnpNull.Rda"))) {
    load(paste0(cache, experiment_name, "dnpNull.Rda"))
    return(null_df)
  }

  #load ppi
  g <- load_ppi(cache = cache, min_score = min_score, ppi = ppi)

  samples <- unique(df$sample_name)

  #grab top n genes per sample to iterate over
  v_df = get_topn(df =df, n_genes = n_genes)

  #get rid of dnp
  df <- df %>% dplyr::select(-.data$dnp)
  #do we want to do a grouped split on sample and then iterate directly through the list? probably
  df_list <- df %>% dplyr::group_by(.data$sample_name) %>%
    dplyr::group_split()

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  null_df <-
    foreach::foreach(j = iterators::iter(df_list),
                     .packages = "disruptr",
                     .combine = "rbind"
    ) %dopar% {


      mem_df <- list()
      agg_df <- list()
      z <- 1
      for(i in 1:n) {

        #get the info on the top genes for a given sample
        j_sample <- j$sample_name[1]
        v_j <- v_df[,j_sample][[1]][[1]] #Indexing is gross but this just grabs the top n genes for a given sample
        #scramble the graph
        g_i <- get_random_graph(g)

        #calc change in network potential
        dnp_ij <- calc_dnp_i(j, g_i, v_rm = v_j, keep_all = FALSE)

        #for testing only
        # dnp_ij <- dplyr::filter(j, .data$gene_name %in% v_j)
        # dnp_ij$dnp <- rnorm(nrow(dnp_ij), mean = 50, sd = 20)

        mem_df[[i]] <- dnp_ij

        if(i %% 10 == 0) {
          agg_df[[z]] <- combine_null(mem_df)
          agg_df[[z]]$z <- z
          mem_df <- list() #clear "mememory
          z <- z + 1#progress ticker
        }
      }
      return(final_combine(agg_df)) #pass back to the main foreach loop

    }
  # null_j_df <- cell_line_df %>%
  #   dplyr::group_by(.data$gene_name, .data$sample_name) %>%
  #   dplyr::summarise(mean_dnp = mean(.data$np),
  #                    sd_dnp = sd(.data$dnp),
  #                    n = dplyr::n())
  #return(null_j_df)



  #close out parallel execution
  parallel::stopCluster(cl)

  #save
  save(null_df, file = paste0(cache, experiment_name, "dnpNull.Rda"))

  #return
  return(null_df)
}

#' Helper function for compute_null - returns a graph with randomly permuted edges.
#'
#' currently just a wrapper for igraph::rewire but may add more functionality in the future
#'
#' @param g graph to be permuted
#'
#' @seealso [igraph::rewire()]
#' @export

get_random_graph <- function(g) {
  g2 <- igraph::rewire(g, igraph::keeping_degseq(niter = igraph::vcount(g)*100))
  return(g2)
}

#' Helper function for compute_null - returns the top n genes by dnp for each sample
#' @inheritParams compute_null
#' @importFrom magrittr %>%

get_topn <- function(df, n_genes) {
  df %>% dplyr::group_by(.data$sample_name) %>%
    dplyr::slice_max(order_by = .data$dnp, n = n_genes) %>%
    dplyr::select(.data$gene_name, .data$sample_name) %>%
    tidyr::pivot_wider(names_from = .data$sample_name,
                       values_from = .data$gene_name,
                       values_fn = list)
}
#' .combine function for compute_null foreach looping structure
#'
#' @param x aggregated data structure
#' @param y task returned by inner foreach loop of compute_null
#'
#' @param ... arbitrary number of tasks returned by the inner foreach loop of
#'            compute_null
#' @importFrom magrittr %>%
#' @export

combine_null <- function(x) {

  #don't want to try and compute sd with one sample... thats where the trouble is ocming I t
  x <- dplyr::bind_rows(x) %>%
    dplyr::group_by(.data$gene_name, .data$sample_name) %>%
    dplyr::summarise(mean_dnp = mean(.data$dnp),
                     sd_dnp = sd(.data$dnp, na.rm = TRUE),
                     n = dplyr::n())

  return(x)
}

#' final .combine function to run in compute_null foreach looping structure
#'
#' @param x aggregated info
#' @param ... arbitrary number of tasks returned by the inner foreach loop of
#'            compute_null
#' @importFrom magrittr %>%
#' @export
#'
#'
final_combine <- function(x) {

  x <- dplyr::bind_rows(x)
  x <- x %>%
    dplyr::group_by(.data$gene_name, .data$sample_name) %>%
    dplyr::summarise(mean_dnp = mean(.data$mean_dnp),
                     sd_dnp = mean(.data$sd_dnp, na.rm = TRUE),
                     n = sum(.data$n))
  return(x)
}
