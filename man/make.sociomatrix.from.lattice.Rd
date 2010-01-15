\name{make.sociomatrix.from.lattice}
\alias{make.sociomatrix.from.lattice}
\title{make.sociomatrix.from.lattice}
\description{
Given a set of points in two dimensions, produce a sociomatrix where a link exists if their distance is exactly 1.
}
\usage{
make.sociomatrix.from.lattice(pts.nby2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pts.nby2}{An n-by-2 matrix whose rows are (x,y) coordinates of points.}
}
\value{
  A sociomatrix whose entries correspond to indicators of whether the points are exactly 1 unit apart.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{
# generate a 5-by-5 grid of points.
pts <- cbind(rep(1:5,5),sort(rep(1:5,5)))
socio <- make.sociomatrix.from.lattice(pts)

}

