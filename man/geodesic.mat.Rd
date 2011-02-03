\name{geodesic.mat}
\alias{shortest.path.maker}
\alias{geodesic.mat}
\title{geodesic.mat}
\description{
Given a k-by-2 or k-by-3 edge list, calculate the minimum distance between any two points using Dijkstra's Algorithm.
}
\usage{
  geodesic.mat(edgelist, node.ids=NULL)
  shortest.path.maker(edgelist, sourcepoint=1, destpoint=NULL, node.ids=NULL)
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
geodesic.mat outputs an n-by-n matrix; cell (i,j) contains the shortest
geodesic distance from point i to point j.

shortest.path.maker outputs a vector of distances from the source node,
accurate as far as the distance for ``destpoint''.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{

latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
lattice.edges <- sociomatrix.to.nby3.edges(make.sociomatrix.from.lattice(latt)$sociomatrix)
latt.geo <- geodesic.mat(lattice.edges)

#Distances from the corner node.
corner.dist <- shortest.path.maker(lattice.edges, 1)

}
