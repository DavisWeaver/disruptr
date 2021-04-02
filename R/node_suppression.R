#' Function to eliminate a node from a network g and calculate the change in some measure of network state
#'
#' this function is still under development.
#'
#' @param g igraph network object
#' @param v_rm index of vertices to remove
#' @param state_function function to use to calculate network state before and after node_repression
#' @param only_neighbors if state function makes calculation for each node, should recalculation be limited to neighbors of the removed node?
#' @param ncores number of cores to use for parallel processing.
#'
#' @param ... additional parameters passed to state function.
#' @export

node_repression <- function(g, v_rm, state_function = calc_np_all, only_neighbors, ncores, ...) {

  if(length(v_rm == 1)) {
    g_new <- igraph::delete_vertices(g = g, v = v_rm)
    s_old <- state_function(g = g, ...)
    s_new <- state_function(g = g_new, ...)
    s_old <- s_old[names(s_new)]
    return(s_new - s_old)
  }
}
