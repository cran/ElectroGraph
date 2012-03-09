\name{dist.coords}
\alias{dist.coords}
\title{dist.coords}
\description{
Solve for the Euclidean distance matrix between two sets of points.
}
\usage{
dist.coords(coords, coords2=NULL)
}
\arguments{
  \item{coords, coords2}{Coordinate matrices of dimension n1-by-k (and
    optionally n2-by-k).}
}
\value{
  If one argument is entered, returns an n1-by-n1 matrix with the
  Euclidean distance between all pairs of points. If two, the matrix is
  the distances between all pairs of points across the two matrices.
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{

coords1 <- array(rnorm(20), c(10,2))
coords2 <- array(rnorm(20), c(10,2))

dist.coords(coords1)
dist.coords(coords1, coords2)

}

