\name{network.components}
\alias{network.components}
\title{network.components}
\description{
Given a sociomatrix, determine the membership of connected components.
}
\usage{
network.components (sociomatrix, pseudo.diameter.bridge=2,
  minimum.relative.strength.for.tie = 1e-08)
}
\arguments{
  \item{sociomatrix}{An n-by-n sociomatrix, where tie strengths indicate inverse distance.}
  \item{pseudo.diameter.bridge}{When constructing a pseudo-diameter
    measurement, this represents the additional distance between
    components.}
  \item{minimum.relative.strength.for.tie}{If a tie value is below this,
  the tie is omitted from the analysis.}
}
\value{
  A list containing the following items:	
  \item{components}{A list whose members are vectors of nodes within each connected component.}
  \item{diameters}{The geodesic diameters of each component.}
  \item{global.pseudo.diameter}{A sum of all diameters in the set plus the psuedo-diameter bridges required to putatively connect all components.}
  \item{component.vector}{The component in which each node is located.}
  \item{geodesic}{The full geodesic path length matrix of the system.}
}
\note{This function is called automatically during the creation of an ElectroGraph object.}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{
sources <- c(1,2,3,4)
sinks <- c(2,3,1,5)
socio <- make.sociomatrix.from.edges(cbind(sources,sinks),symmetric=TRUE)

comps <- network.components(socio)

}

