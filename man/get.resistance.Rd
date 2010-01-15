\name{get.resistance}
\alias{get.resistance}
\alias{compute.volts.asymmetric}
\alias{compute.volts.symmetric}
\title{Calculating electro-social properties of a network}
\description{
Given a pair of nodes to act as source and sink, determine the resulting current flow.
}
\usage{
compute.volts.symmetric (sociomatrix,src,snk)
compute.volts.asymmetric (sociomatrix,src,snk)
get.resistance (e.graph, sources, sinks)
}
\arguments{
  \item{sociomatrix}{An n-by-n sociomatrix.}
  \item{e.graph}{An electrograph object.}
  \item{src, snk}{Single values for the source and sink.}
  \item{sources, sinks}{Vectors of source and sink values.}

}
\value{
	Returns a list containing the following items: 
	\item{resist.eq}{The equivalent resistance between the source and sink.}
	\item{voltages}{Each column contains the electric potential at each node given a source/sink pair.}
	\item{current.out}{The total current flowing into and out of each node given a source/sink pair.}
	\item{reference.mean}{As a check, computed mean of the sociomatrix.}
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\note{
	The ``compute'' functions are subsumed within get.resistance() and are given for diagnostic purposes. In general, the full execution of ``get.resistance'' is more reliable and faster; this is called in general form by ``electrograph.exam''.
}
\examples{

latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
latt.mat <- make.sociomatrix.from.lattice(latt)
print(compute.volts.symmetric(latt.mat,1,25))
latt.mat[24,25] <- 0
print(compute.volts.asymmetric(latt.mat,1,25))
print(compute.volts.asymmetric(latt.mat,25,1))

latt.e <- electrograph (latt.mat)
print(get.resistance(latt.e,c(1,2,3),c(23,24,25)))


}

