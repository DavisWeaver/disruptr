% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_null.R
\name{get_random_graph}
\alias{get_random_graph}
\title{Helper function for compute_null - returns a graph with randomly permuted edges.}
\usage{
get_random_graph(g)
}
\arguments{
\item{g}{graph to be permuted}
}
\description{
currently just a wrapper for igraph::rewire but may add more functionality in the future
}
\seealso{
\code{\link[igraph:rewire]{igraph::rewire()}}
}
