\name{ElectroGraph-package}
\alias{ElectroGraph-package}
\alias{ElectroGraph}
\docType{package}
\title{
The ElectroGraph package: Analyzing and visualizing relational data
}
\description{
ElectroGraph takes as input an edge list or sociomatrix, either valued or binary, and computes various statistics for node and edge importance.
}
\details{
\tabular{ll}{
Package: \tab ElectroGraph\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-05-15\cr
License: \tab GPL, v2\cr
LazyLoad: \tab yes\cr
}
To create an ElectroGraph object, specify an edge list or sociomatrix.
}
\author{
Andrew C. Thomas

Maintainer: Andrew C. Thomas <acthomas@fas.harvard.edu>
}
\references{
Andrew C. Thomas (2009) ``Hierarchical Models for Relational Data''.
}
\keyword{ package }
\examples{
latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
latt.e <- electrograph (make.sociomatrix.from.lattice(latt))
plot(latt.e)
}
