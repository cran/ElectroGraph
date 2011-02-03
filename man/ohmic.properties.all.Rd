\name{ohmic.properties.all}
\alias{ohmic.properties.all}
\title{ohmic.properties.all}
\description{
Compute a full analysis of Ohmic properties of the system.
}
\usage{
ohmic.properties.all(e.graph, sample.edges=FALSE, sample.per=NULL, verbose=FALSE)
}
\arguments{
  \item{e.graph}{An object initialized by the electrograph() function.}
  \item{sample.edges}{If TRUE, samples of edges are taken as opposed to the full set.}
  \item{sample.per}{Sets the number of samples taken if sample.edges is TRUE. Defaults to the square root of the number of edges.}
  \item{verbose}{Enable additional (diagnostic) output.}
}

\value{
The original electrograph object, augmented with electro-social quantities:

  \item{ohmic.distance.mat}{The equivalent electro-social resistances for each pair of nodes in matrix form; if sample.edges is TRUE, edge importance is approximated.}
  \item{conductances}{A vector of electro-social conductances for each specified pair of nodes.}
  \item{voltages}{A matrix with each column representing the voltages at each node.}
  \item{currents.node}{A matrix with each column representing the current through each node.}
  \item{source.sink}{The source-sink vectors used in the examination.}
  \item{current.matrix}{The average current through each directed edge.}
}
\note{
  This is currently performed automatically when the constructor function ``electrograph'' is called, but can be executed manually as well.
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{
latt <- cbind(rep(1:5,5),sort(rep(1:5,5)))
latt.e <- electrograph (make.sociomatrix.from.lattice(latt)$sociomatrix, ohmic=FALSE)
latt.exam <- ohmic.properties.all(latt.e)
}
