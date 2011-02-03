\name{get.resistance}
\alias{get.resistance}
\alias{compute.volts.asymmetric.socio}
\alias{compute.volts.asymmetric.edges}
\alias{compute.volts.symmetric.socio}
\alias{compute.volts.symmetric.edges}
\title{Calculating Ohmic properties of a network}
\description{
Given a pair of nodes to act as source and sink, determine the resulting current flow.
}
\usage{
compute.volts.symmetric.edges (edgelist, src, snk, fidelities = NULL)
compute.volts.asymmetric.edges (edgelist, src, snk, fidelities = NULL)
compute.volts.symmetric.socio (sociomatrix, src, snk, fidelities = NULL)
compute.volts.asymmetric.socio (sociomatrix, src, snk, fidelities = NULL)
get.resistance (e.graph, sources, sinks, verbose=FALSE)
}
\arguments{
  \item{sociomatrix}{An n-by-n sociomatrix.}
  \item{edgelist}{A k-by-\{2,3,4\} edge list.}
  \item{e.graph}{An electrograph object.}
  \item{src, snk}{Single values for the source and sink.}
  \item{sources, sinks}{Vectors of source and sink values.}
  \item{fidelities}{If set, allows for ties that are perceived to be
    less than ``friends''. (Use with caution -- results have not been peer reviewed!}
  \item{verbose}{Enable additional (diagnostic) output.}

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
latt.mat <- make.sociomatrix.from.lattice(latt)$sociomatrix

print(compute.volts.symmetric.socio(latt.mat,1,25))
latt.mat[24,25] <- 0
print(compute.volts.asymmetric.socio(latt.mat,1,25))
print(compute.volts.asymmetric.socio(latt.mat,25,1))

latt.e <- electrograph (latt.mat)
print(get.resistance(latt.e,c(1,2,3),c(23,24,25)))


}

