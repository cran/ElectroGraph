\name{geodesic.mat}
\alias{geodesic.mat}
\title{geodesic.mat}
\description{
Given a sociomatrix corresponding to social distance, calculate the minimum distance between any two points using the Floyd-Warshall Algorithm.
}
\usage{
geodesic.mat(sociomatrix)
}
\arguments{
  \item{sociomatrix}{An n-by-n sociomatrix corresponding to tie strength (or inverse distance).}
}
\value{
An n-by-n matrix; cell (i,j) contains the shortest distance from point i to point j.
}
\note{
This function is automatically executed by the electrograph() initialization routine.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{

latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
latt.mat <- make.sociomatrix.from.lattice(latt)$sociomatrix
latt.geo <- geodesic.mat(latt.mat)

}
