\name{network.projection}
\alias{network.projection}
\title{network.projection}
\description{
Given a table of ideal distances and a connectivity matrix, calculate the projected coordinates under force-directed placement.
}
\usage{
network.projection(connectivity,
                   force.mode=c("fruchterman.reingold","kamada.kawai"),
                   layout.dimension=2, 
		   ego.focus=NULL, 
		   initial.coordinates = NULL,
                   verbose=1,
                   projection.iterations=5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{connectivity}{An n-by-n matrix containing the values of the
    connective strengths among the n individuals (inverse distance).}
  \item{force.mode}{Either of the above two standard methods for force-directed placement.}
  \item{layout.dimension}{Of 1, 2 or 3, the dimension of the output coordinates.}
  \item{ego.focus}{A vector of node indices for which priority should be made highest for their relative distances.}
  \item{initial.coordinates}{If specified, the algorithm will begin at the specified coordinates rather than a randomly generated set.}
  \item{verbose}{If TRUE, display additional information while running.}
  \item{projection.iterations}{The maximum number of cycles that the
    projection algorithm will run.}
	
}
\value{
	An n-by-(layout.dimension) table of coordinates.
}
\note{
	The algorithm can technically handle any dimensionality in the data, and is written for expandability. However, in its current form, as only two-dimensional plotting is supported in the package (and three-dimensional visualization the maximum supported by the human brain) the layout dimension is restricted in this routine.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{

latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
latt.e <- electrograph (make.sociomatrix.from.lattice(latt)$sociomatrix,
shortest.paths=TRUE)
latt.coords <- network.projection(1/latt.e$geodesic)

print(latt.coords)

}

