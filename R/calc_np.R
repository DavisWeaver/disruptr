#' function to calculate the network potential for each protein in a user-provided vector
#'
#' @param exp expression vector - assumed to be a named vector where the values are expression and the names are the gene name
#' @param g igraph object - will be filtered so that only nodes found in both exp and g are kept
#' @param v character vector of nodes over which to calculate network potential.
#' @param neighbors named list containing the neighbors for each node of graph g. If not provided,
#'        it will be computed
#' @return dataframe containing network potential for each of the inputed gene names.
#'

calc_np_all_legacy <- function(exp, g, v = as.character(names(igraph::V(g))), neighbors = NULL) {

  #first add expression to the subgraph.
  g <- add_expression(exp = exp, g = g)

  vertices <- as.character(names(igraph::V(g))) #in most cases this will be the same as `v`

  # remove any names of exp that are not in the graph.
  exp <- exp[names(exp) %in% vertices]

  # need to do the same thing for v (because sometimes there will be neighbors of a given node that aren't in the ppi)
  v <- v[v %in% vertices]

  #get a list of neighbors for each node
  if(is.null(neighbors)) {
    neighbors <-
      lapply(v,
             get_neighbors, g = g)
    names(neighbors) <- v

  }
  #loop over all vertices
  np_vec <- vector(mode = "numeric", length = length(v))
  names(np_vec) <- v
  vertex_list <- igraph::V(g) #slightly different than "vertices" - used for indexing
  for(i in v) {
    neighbors_named <- as.character(names(vertex_list[neighbors[[i]]])) #grab named vec of neighbors for each vertex
    c_j <- sum(exp[neighbors_named]) #sum up the concentration of all neighbors
    c_i <- exp[i]
    np_vec[i] <- calc_np(c_i = c_i, c_j = c_j)
  }

  #gonna do some data wrangling to get the vertices in the same order as the input expression vector
  if(length(np_vec) < length(exp)) {
    exp <- exp[names(exp) %in% names(np_vec)]
  }
  np_vec <-np_vec[names(exp)]

  return(np_vec)
}

#' function to calculate the network potential for each protein in a user-provided vector - cpp internal version
#'
#' @param exp expression vector - assumed to be a named vector where the values are expression and the names are the gene name
#' @param g igraph object - will be filtered so that only nodes found in both exp and g are kept
#' @param v character vector of nodes over which to calculate network potential.
#' @param neighbors named list containing the neighbors for each node of graph g. If not provided,
#'        it will be computed
#' @return dataframe containing network potential for each of the inputed gene names.
#'
#' @export

calc_np_all <- function(exp, g, v = as.character(names(igraph::V(g))),
                         neighbors = NULL) {

  #first add expression to the subgraph.
  g <- add_expression(exp = exp, g = g)

  vertices <- as.character(names(igraph::V(g))) #in most cases this will be the same as `v`

  # remove any names of exp that are not in the graph.
  exp <- exp[names(exp) %in% vertices]

  # need to do the same thing for v (because sometimes there will be neighbors of a given node that aren't in the ppi)
  # actually how is that even possible? leaving it in for now
  v <- v[v %in% vertices]
  #v we will just use as a logical flag

  #re-order exp to have the same order as vertices
  exp <- exp[vertices]

  if(is.null(neighbors)) {
    neighbors <-
      lapply(vertices,
             get_neighbors, g = g)
    names(neighbors) <- vertices
  } else {
    neighbors <- neighbors[vertices] #index to make sure we don't have any list items we don't need and to make sure they are in the same order as exp and vertices
  }

  #run cpp function to do the actual calculation on each node
  np_vec <- fcalc_np_all(neighbors = neighbors, vertices = vertices, v = v, exp = exp)
  #gonna do some data wrangling to get the vertices in the same order as the input expression vector
  # if(length(np_vec) < length(exp)) {
  #   exp <- exp[names(exp) %in% names(np_vec)]
  # }

  return(np_vec)
}

#' calculate network potential for one node.
#'
#' @param c_i expression for a given node.
#' @param c_j vector of expressions for each neighbor of c_i
#'
#' @export

calc_np <- function(c_i, c_j) {
  g_i <- c_i*log(c_i/(sum(c_j)))
  return(g_i)
}


#' attach expression values from user-provided expression vector to graph.
#'
#' @inheritParams calc_np_all
#'
#' @return subgraph of g containing only shared keys with exp and with expression attached.

add_expression <- function(exp, g) {
  #subset exp and g so they contain the same index values
  vertices <- as.character(names(igraph::V(g)))
  keep_vertices <- vertices[vertices %in% names(exp)]
  exp <- exp[keep_vertices]

  #create new subgraph
  g <- igraph::induced_subgraph(g, keep_vertices)

  #attach expression as an attribute of that subgraph
  g <- igraph::set_vertex_attr(g, name = "expression",value = exp)

  return(g)
}

#' function to get graph neighbors (along with their expression values) for a given gene in a given network g
#'
#' just a wrapper around [igraph::neighbors()] for convenience
#'
#' @param gene gene to grab neighbors from.
#' @inheritParams calc_np_all
#'
#' @return named numeric vector.

get_neighbors <- function(gene, g) {
  neighborGenes <- as.vector(igraph::neighbors(g, gene))
  return(unique(neighborGenes))
}
