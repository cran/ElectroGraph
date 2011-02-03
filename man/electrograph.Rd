\name{electrograph}
\alias{electrograph}
\title{electrograph}
\description{
Given a sociomatrix or an edge list, produce an electrograph object for display and analysis.
}
\usage{
electrograph(input, make.edgelist.symmetric=TRUE,
             solve.for.shortest.paths=TRUE, solve.for.betweenness=FALSE,
             solve.for.canadian.betweenness=FALSE,
	     ohmic.properties=FALSE, fidelities=NULL, 
	     verbose=FALSE, substitute.names=NULL, ...)
}
\arguments{
  \item{input}{A sociomatrix (n-by-n) or edge list (k-by-(2,3,4).}
  \item{make.edgelist.symmetric}{Whether the input edge list should be taken as symmetric.}
  \item{solve.for.shortest.paths}{If TRUE, calculates the geodesic
    distances between nodes.}
  \item{solve.for.betweenness}{If TRUE, calculates the standard geodesic
    betweenness centrality for all edges.}
  \item{solve.for.canadian.betweenness}{If TRUE, calculates the Canadian
    betweenness centrality for all edges.}
  
  \item{ohmic.properties}{If TRUE, calculates the Ohmic properties of the network nodes and ties through electrograph.exam. This will be lengthy if the system has a large number of nodes.}
  \item{fidelities}{If set, allows for ties that are perceived to be
    less than ``friends''. (Use with caution -- results have not been peer reviewed!}
  \item{verbose}{Additional stage-by-stage feedback is given if TRUE.}
  \item{substitute.names}{A k-by-2 matrix with one row containing original node designations and the other containing replacement names}
  \item{...}{Additional options to be passed to electrograph.exam.}
}
\details{
  The input object can have one of four forms:
    -A square sociomatrix.
    -An n-by-2 matrix containing social arcs; the individuals in column 1 consider those in column 2 to be friends with social strength 1. By default, this is symmetric.
    -An n-by-3 matrix containing social arcs; the individuals in column 1 consider those in column 2 to be friends with social strength as described in column 3. By default, this is symmetric.
    -An n-by-4 matrix containing dyadic arcs. The connection (col1,col2) is described in column 3; the reciprocal connection (col2,col1) is described in column 4. This by nature overrides the option ``symmetric'' to be false unless the input is truly symmetric.
}
\value{
  A list, of class ``electrograph'' containing the following items:	
  \item{grand.sociomatrix}{The full sociomatrix of the system.}
  \item{grand.sociomatrices}{A list of sociomatrices of the components of the system.}
  \item{component.vector}{The component in which each node is located.}
  \item{geodesic}{The full geodesic path length matrix of the system.}
  \item{diameters}{The geodesic diameters of each component.}
  \item{global.pseudo.diameter}{A sum of all diameters in the set plus the psuedo-diameter bridges required to putatively connect all components.}
  \item{symmetric}{Is the resulting sociomatrix symmetric?}
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\source{
Andrew C. Thomas (2010). ``Ohmic Circuit Interpretations of Network
  Distance and Centrality''. Unpublished: available at
  http://www.acthomas.ca/.
Andrew C. Thomas and Stephen E. Fienberg (2011). ``Exploring the Consequences of IED Deployment with a Generalized Linear Model Implementation of the Canadian Traveller Problem''. http://arxiv.org/abs/1012.4185
}
\examples{

sources <- c(1,2,3,4)
sinks <- c(2,3,1,5)
socio <- electrograph(cbind(sources,sinks))

}

