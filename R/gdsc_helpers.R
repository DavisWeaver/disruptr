#' function to calculate network potential for every gene/cell line combo in the GDSC
#'
#' @param cache
#'
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#'
#' @export

calc_gdsc_np <- function(cache = NULL, ncores = 1) {
  if(is.null(cache)) {
    stop("Please specify a filepath specifying where to store/retrieve data")
  }
  #just return the processed data if its already been done
  if(file.exists(paste0(cache, "/gdsc_np.Rda"))) {
    load(paste0(cache, "/gdsc_np.Rda"))
    message("using cached processed data")
    return(df_np)
  }
  #grab ppi to use for the scaffold
  g <- crosstalkr::prep_stringdb(cache = cache, min_score = 400)

  #clean gdsc
  df_list <- clean_gdsc(cache = cache)

  #grab expression data
  df_exp <- df_list[[3]]

  #get list of cell lines to iterate over
  cell_lines <- unique(df_exp$cell_line)

  #set up parallel execution
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  out_list <-
    foreach::foreach(i = 1:length(cell_lines)) %dopar% {
      #isolate one cell line
      df_i <- df_exp %>%
        dplyr::filter(cell_line == cell_lines[i])

      #add np to df_i for that cell line
      disruptr::calc_np_i(df_i, g = g)

    }
  #close out parallel execution
  parallel::stopCluster(cl)
  df_np = dplyr::bind_rows(out_list)

  save(df_np, file = paste0(cache, "/gdsc_np.Rda"))
  return(df_np)
}

#' helper function to calculate np for one cell line in the gdsc
#'
#' @param df dataframe with one cell line + log expression
#' @param g igraph object containing ppi info
#'
#' @return same dataframe with np calculated for each gene.
#'
#' @export

calc_np_i <- function(df, g) {
  #grab expression vector
  exp <- df$log_expression
  names(exp) <- df$gene_symbols

  #calc np for that expression vector
  np <- calc_np_all(exp, g)
  np_df <- data.frame(gene_symbols = names(np),
                      np = np)


  df <- dplyr::left_join(df, np_df)
  return(df)
}

#' function to do in-silico repression for a given gene on all cell lines in the GDSC
#'
#' @inheritParams calc_gdsc_np
#' @param v_rm which node should we remove
#'
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#'
#' @export


GDSC_repress <- function(cache = NULL, v_rm = "EGFR", ncores = 1) {

  if(is.null(cache)) {
    stop("Please specify a filepath specifying where to store/retrieve data")
  }

  if(file.exists(paste0(cache, "/", v_rm, "_gdsc_np.Rda"))) {
    load(paste0(cache, "/", v_rm, "_gdsc_np.Rda"))
    message("using cached processed data")
    return(df_np)
  }

  #get PPI
  g <- crosstalkr::prep_stringdb(cache = cache, min_score = 400)

  #get or compute network potential for gdsc
  df <- calc_gdsc_np(cache = cache, ncores = ncores)

  #compute the change in np for EGFR for each cell line.
  cell_lines <- unique(df$cell_line)

  #set up parallel computing
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  out_list <-
    foreach::foreach(i = 1:length(cell_lines),
                     .packages = "disruptr") %dopar% {

                       #get named expression vector
                       exp <- calc_exp(df = df, cell_lines = cell_lines, i = i)

                       #calculate change in np when v_rm is removed
                       repression_vector <- node_repression(g = g, v_rm = v_rm,
                                                            state_function = calc_np_all,
                                                            exp = exp,
                                                            neighbors_only = TRUE)

                       #clean up and bind
                       bind_repression(df = df, cell_lines = cell_lines, i = i,
                                       v_rm = v_rm,
                                       repression_vector = repression_vector)
                     }
  parallel::stopCluster(cl)

  #bind the data frame back together
  df_np = dplyr::bind_rows(out_list)

  #save the output
  save(df_np, file = paste0(cache, "/", v_rm, "_gdsc_np.Rda"))

  #exit function
  return(df_np)

}

#' helper function to grab the named expression vector for a given cell line
#'
#' @param df tidy dataframe containing expression data for a number of cell lines
#' @param cell_lines unique vector of cell lines
#' @param i index value for `cell_lines
#'
#' @export

calc_exp <- function(df, cell_lines, i) {
  df_i <- dplyr::filter(df, cell_line == cell_lines[i])
  exp <- df_i$log_expression
  names(exp) <- df_i$gene_symbols
  return(exp)
}

#' helper function to bind the repression vector to df_i in gdsc_repress
#'
#' @inheritParams calc_exp
#' @param repression_vector
#'
#' @export

bind_repression <- function(df, cell_lines, i, repression_vector, v_rm) {
  df_i <- dplyr::filter(df, cell_line == cell_lines[i])
  df_repression <- data.frame(gene_symbols = names(repression_vector),
                              np_diff = repression_vector)

  df_i <- dplyr::left_join(df_i, df_repression)
  df_i <- dplyr::mutate(df_i,
                        np_diff = ifelse(gene_symbols == v_rm, abs(np), np_diff),
                        np_diff = ifelse(is.na(np_diff), 0, np_diff))
  return(df_i)

}

#' function to summarise output from `GDSC_repress`
#'
#' @param df data frame containing NP difference from in-silico repression of a given gene
#'
#' Computes the total network potential/ difference in network potential with repression of a given gene,
#' expression for the repressed gene, network potential for the repressed gene.
#'
#' @importFrom magrittr %>%
#'
#' @export

clean_output <- function(df) {
  df <- df %>%
    dplyr::filter(!is.na(np),
                  !is.infinite(np)) %>%
    dplyr::group_by(cell_line) %>%
    dplyr::summarise(total_np = sum(.data$np),
                     diff_np = sum(.data$np_diff),
                     egfr_exp = sum(ifelse(
                       gene_symbols == 'EGFR',
                       .data$log_expression,
                       0)),
                     egfr_np = sum(ifelse(
                       gene_symbols == 'EGFR',
                       .data$np,
                       0))) %>%
    dplyr::mutate(diff_np_scaled = .data$diff_np/.data$total_np)
}



