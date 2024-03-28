####networkstat.Version1.0.0###
###Author:Wang Ningqi####
####度数量 total degree###
#' Igraph network statistics
#'
#' @param input Igraph graph object
#'
#' @return network statistics
#' @export
#' @author  Wang Ningqi <2434066068@qq.com>
#' @examples network_stat(igraph).See details in \code{\link{network_analysis}}
network_stat=function(input){require(igraph)
  cat("Total degree:",sum(igraph::degree(input))) #totaldegree=sum(degree(igraph))#
  #  The size of the graph (number of edges)
  cat("
      Total edges/links:",length(E(input)))##num.edges = length(E(igraph1)) ##
  #  Order (number of vertices) of a graph
  cat("
Total vertices:",length(V(input)))#num.vertices = length(V(igraph1))#
  # connectance
  cat("
Connectance:",edge_density(input,loops=FALSE))
  # (Average degree) degree/vertices
  cat("
Average degree:",mean(igraph::degree(input)))
  # Diameter
  cat("
Diameter:",diameter(input, directed = FALSE, unconnected = TRUE, weights = NULL))
  # Average path length
  cat("
Average path length:",average.path.length(input))
  # Clustering coefficient
  cat("
Global clustering coefficient:",transitivity(input))
  cat("
Number of clusters:",no.clusters(input))
  # Betweenness centralization
  cat("
Betweenness centralization:",centralization.betweenness(input)$centralization)
  # Degree centralization
  cat("
Degree centralization:",centralization.degree(input)$centralization,"\n")}
