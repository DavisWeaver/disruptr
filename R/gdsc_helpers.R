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
