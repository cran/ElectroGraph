\name{network.layout}
\alias{network.layout}
\title{network.layout}
\description{
Given a table of ideal distances and a sociomatrix, calculate the projected coordinates under force-directed placement.
}
\usage{
network.layout(sociomatrix, distances,
               plot.mode=c("fruchterman.reingold","kamada.kawai"),
               layout.dimension=2, ego.focus=NULL, initial.coordinates = NULL,
               niter = 500, initemp = 10, sigma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sociomatrix}{An n-by-n matrix containing the values of ties between n individuals.}
  \item{distances}{An n-by-n matrix containing the ideal distances (typically, by shortest path or electro-social) between n individuals.}
  \item{plot.mode}{Either of the above two methods for force-directed placement.}
  \item{layout.dimension}{Of 1, 2 or 3, the dimension of the output coordinates.}
  \item{ego.focus}{A vector of node indices for which priority should be made highest for their relative distances.}
  \item{initial.coordinates}{If specified, the algorithm will begin at the specified coordinates rather than a randomly generated set.}
  \item{niter, initemp, sigma}{Quantities for use in the simulated annealing portions of the algorithm.}
}
\value{
	An n-by-(layout.dimension) table of coordinates.
}
\note{
	The algorithm can technically handle any dimensionality in the data, and is written for expandability; however, in its current form, as only two-dimensional plotting is supported in the package (and three-dimensional visualization the maximum supported by the human brain) the layout dimension is restricted in this routine.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{

latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
latt.e <- electrograph (make.sociomatrix.from.lattice(latt))
latt.coords <- network.layout(latt.e$grand.socio, latt.e$geodesic)

print(latt.coords)

}

