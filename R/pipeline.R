####### Functions to help with users trying to include these network analyses in a pipeline ######
# These functions will allow users to provide an expression matrix / tidy data frame and then return a tidy data frame with network potential/ change in np

#' main function to compute np from a user-provided expression matrix.
#'
#' @param cache user-provided filepath for where to store data etc
#'
#'
#' @return tidy data frame with one column for expression and another for np
#' @importFrom magrittr %>%
#' @export
#'
#'

compute_np <- function(cache = NULL) {

}

#' helper function to convert expression matrix to tidy dataframe (if not already)
#'
#' @inheritParams compute_np
#'
#' @importFrom magrittr %>%

tidy_expression <- function() {

}
