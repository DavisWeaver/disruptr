#' Function to eliminate a node from a network g and calculate the change in some measure of network state
#'
#' this function is still under development.
#'
#' @param g igraph network object
#' @param v_rm index of vertices to remove
#' @param state_function function to use to calculate network state before and after node_repression
#' @param neighbors_only logical designating whether state function should be calculated for all nodes or just neighbors
#' @param ... additional parameters passed to state function.
#' @export

node_repression <- function(g, v_rm, state_function = calc_np_all,
                            neighbors_only = TRUE, ...) {

  if(length(v_rm == 1)) {
    #define v to be either all nodes or neighbors of v_rm
    if(neighbors_only == TRUE) {
      v <- get_neighbors(v_rm, g) #get neighbors
      v <- as.character(names(igraph::V(g)[v]))
    } else {
      v = as.character(names(igraph::V(g)))
    }
    g_new <- igraph::delete_vertices(g = g, v = v_rm)
    s_old <- state_function(g = g, ...)
    s_new <- state_function(g = g_new, v = v, ...)
    s_old <- s_old[names(s_new)]
    return(s_new - s_old)
  }
}
