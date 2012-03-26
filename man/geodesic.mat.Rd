\name{geodesic.mat}
\alias{geodesic.selected.from.edgelist}
\alias{geodesic.matrix.from.edgelist}
\title{geodesic.mat}
\description{
Given a k-by-2 or k-by-3 edge list, calculate the minimum geodesic distance between any two points using Dijkstra's Algorithm.
}
\usage{
  geodesic.matrix.from.edgelist(edgelist, node.ids=NULL)
  geodesic.selected.from.edgelist(edgelist, sourcepoint=1, destpoint=NULL, node.ids=NULL)
}
\arguments{
  \item{edgelist}{A k-by-2 or k-by-3 matrix containing edge
    information. The first two columns represent the nodes in the arc;
    if included, the third refers to the strength of the tie (inverse
    distance).}
  \item{sourcepoint}{The number of the source node.}
  \item{destpoint}{The farthest node whose distance from the source should
    be calculated. }
  \item{node.ids}{Labels for the nodes in the system.}
}
\value{
geodesic.matrix.from.edgelist outputs an n-by-n matrix; cell (i,j) contains the shortest
geodesic distance from point i to point j.

geodesic.selected.from.edgelist outputs a vector of distances from the source node,
accurate as far as the distance for ``destpoint''.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{

latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
lattice.edges <- sociomatrix.to.nby3.edges(make.sociomatrix.from.lattice(latt)$sociomatrix)
latt.geo <- geodesic.matrix.from.edgelist(lattice.edges)

#Distances from the corner node.
corner.dist <- geodesic.selected.from.edgelist(lattice.edges, 1)

}
