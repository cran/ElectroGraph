\name{betweenness.centralities}
\alias{single.source.betweenness}
\alias{betweenness.centralities}
\alias{canadian.betweenness.one}
\alias{canadian.betweenness.source}
\alias{canadian.betweenness.full}
\title{Calculating Geodesic Betweenness Properties of a Network}
\description{
Functions that
}
\usage{
single.source.betweenness (edgelist, sourcepoint=1, destpoint=NULL, node.ids=NULL)
betweenness.centralities (edgelist,
  path.weight=c("constant","closeness"), verbose=TRUE, node.ids=NULL)

canadian.betweenness.one (edgelist, sourcepoint=1, destpoint=NULL, penalty=20, node.ids=NULL) 
canadian.betweenness.source (edgelist, sourcepoint=1, penalty=20, node.ids=NULL) 
canadian.betweenness.full (edgelist, penalty=20, path.weight=c("constant","closeness"), node.ids=NULL) 

}
\arguments{
  \item{edgelist}{A k-by-2 or k-by-3 matrix containing edge
    information. The first two columns represent the nodes in the arc;
    if included, the third refers to the strength of the tie (inverse
    distance).}
  \item{sourcepoint}{The number of the source node.}
  \item{destpoint}{The number of the destination node.}
  \item{verbose}{Display additional output if TRUE.}
  \item{penalty}{The additional cost to be considered in Canadian Betweenness if
    the route is blocked.}
  \item{path.weight}{Should the betweennesses of edges be calculated with all
    source-sink pairs equal, or weighted inversely by their geodesic
    distance?}
  \item{node.ids}{Labels for the nodes in the system.}
}
\value{
  Objects output the betweenness measures for the edges in the system,
  given the sources and destinations as inputs.
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{

  #Ring lattice.
  edgelist <- cbind(c(1:5, 2:5, 1), c(2:5,1, 1:5))

  shortest.path.maker (edgelist)
  single.source.betweenness (edgelist)
  betweenness.centralities (edgelist, path.weight="constant")

  canadian.betweenness.full (edgelist)
}

