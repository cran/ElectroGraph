\name{summary.electrograph}
\alias{summary.electrograph}
\alias{print.electrograph}
\title{summary.electrograph}
\description{
Given an electrograph object, summarize the qualities of its nodes.
}
\usage{
\method{summary}{electrograph}(object, ...)
\method{print}{electrograph}(x, ...)
}
\arguments{
  \item{x, object}{An object of class ``electrograph''.}
  \item{...}{(no further options currently available)}
}
\value{
	"print" returns an abbreviated list of contents of the object.
	"summary" Returns a list containing several measures:
  \item{in.deg, out.deg}{The sums of the inbound and outbound weights of each tie from each node.}
  \item{cent.geodesic}{Geodesic centrality: the sums of the inverses of the minimum distances to each other node in the system.}
  \item{cent.electro}{Electro-social centrality: the sums of the effective conductances from a node to each other node in the system. Defined only if electrograph.exam has been run on the object.}
  \item{current.cent}{Electro-betweenness measure: the average current through each node, given each node pair for which the effective conductance was measured.}
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{

sources <- c(1,2,3,4)
sinks <- c(2,3,1,5)
socio <- electrograph(cbind(sources,sinks))

print(socio)
summary(socio)

}



