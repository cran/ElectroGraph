\name{sociomatrix.to.nby3.edges}
\alias{sociomatrix.to.nby3.edges}
\title{sociomatrix.to.nby3.edges}
\description{
Given a sociomatrix, produce an edge/arc list whose columns are the source,
destination and value of each arc.
}
\usage{
sociomatrix.to.nby3.edges(sociomatrix, keep.labels=FALSE)
}
\arguments{
  \item{sociomatrix}{An n-by-n matrix whose elements are the strengths
    of arcs corresponding to the row and column indices.}
  \item{keep.labels}{If true, maintains the row names as the node
    labels; if false, assigns node IDs as 1:nrow(sociomatrix).}
}
\value{
A k-by-3 matrix, whose columns are the source, destination and value of
each of k arcs.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{

latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
lattice.socio <- make.sociomatrix.from.lattice(latt)$sociomatrix
lattice.edges <- sociomatrix.to.nby3.edges(lattice.socio)

}
