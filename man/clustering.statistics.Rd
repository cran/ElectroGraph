\name{clustering.statistics}
\alias{clustering.statistics}
\title{clustering.statistics}
\description{
Calculate the transitive counts and cycles in a sociomatrix.
}
\usage{
   clustering.statistics(sociomatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sociomatrix}{An n-by-n sociomatrix.}
}
\value{
  A data frame containing the fraction of possible cycles and transitive triples
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{
# generate a 5-by-5 grid of points.
pts <- cbind(rep(1:5,5),sort(rep(1:5,5)))
socio <- make.sociomatrix.from.lattice(pts)
clustering.statistics(socio$sociomatrix)

}

