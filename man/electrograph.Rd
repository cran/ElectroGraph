\name{electrograph}
\alias{electrograph}
\title{electrograph}
\description{
Given a sociomatrix or an edge list, produce an electrograph object for display and analysis.
}
\usage{
electrograph(input, symmetric=TRUE, perform.exam=TRUE,
					enemies.allowed=FALSE, ...)
}
\arguments{
  \item{input}{A sociomatrix (n-by-n) or edge list (k-by-(2,3,4).}
  \item{symmetric}{Whether the input edge list should be taken as symmetric.}
  \item{perform.exam}{If true, automatically conducts a full electrograph.exam of the created object. 
This will be lengthy if the system has a large number of nodes.}
  \item{enemies.allowed}{If enabled, negative values can be considered as tie strengths. (Use with caution -- results have not been peer reviewed!}
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
\author{Andrew C. Thomas <acthomas@fas.harvard.edu>}
\examples{

sources <- c(1,2,3,4)
sinks <- c(2,3,1,5)
socio <- electrograph(cbind(sources,sinks))

}

