\name{plot.electrograph}
\alias{plot.electrograph}
\title{plot.electrograph}
\description{
Given an electrograph object, display a two-dimensional representation.
}
\usage{
\method{plot}{electrograph}(x, distance.mode = c("shortest.path","electro.social"), 
					plot.mode = c("kamada.kawai", "fruchterman.reingold"), 
					ego.focus = NULL, manual.coords = NULL, redo.coordinates = FALSE,
					just.coordinates = FALSE, node.colors = NULL,
					label.colors = NULL, pts.cex = 1.5, label.cex = 1.5,
					component.border.col = 8, 
					edge.thickness = c("standard", "electro.betweenness"), 
					max.thick = 2,	source.sink.pair = NULL,
					previous.electrograph.plot.object = NULL, 
					tick.marks=FALSE,
					bound.size=1,
					width.height.factor=1,
						...)
}
\arguments{
  \item{x}{An object of class ``electrograph''.}
  \item{distance.mode}{Should the projected distances be calculated using the shortest path, or the equivalent conductance of the electro-social approach?}
  \item{plot.mode}{Which force-directed placement method should be used to find node coordinates?}
  \item{ego.focus}{A vector of nodes whose importance should be considered highest in the plot.}
  \item{manual.coords}{If desired, the manual input of coordinates for each node.}
  \item{redo.coordinates}{If true, recalculate the true distances between points using the desired distance method.}
  \item{just.coordinates}{If true, do not output a plot.}
  \item{node.colors}{The colors of the nodes to be plotted.}
  \item{label.colors}{The colors of node labels.}
  \item{pts.cex}{Point sizes.}
  \item{label.cex}{Label sizes.}
  \item{component.border.col}{The color assigned to the boxes separating graph components.}
  \item{edge.thickness}{Determines whether the sociomatrix or electro-social betweenness will determine the thickness of edges in the plot.}
  \item{max.thick}{The maximum thickness of an edge in plotting.}
  \item{source.sink.pair}{A vector of length 2 that specifies an electro-social current path to be highlighted.}
  \item{previous.electrograph.plot.object}{A previously produced graph object, for the sake of lining up a new plot with a previous one.}
  \item{tick.marks}{If TRUE, place tick marks at equal distances on each axis.}
  \item{bound.size}{Denotes the degree extra space to surround each component in the plot.}
  \item{width.height.factor}{The desired ratio of width to height in the resulting plot.}
  \item{...}{Additional options to pass to ``plot''.}
  
}
\value{
	Returns an object of class ``electrograph.plot'' which can be replotted or used in later analysis.
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{

sources <- c(1,2,3,4)
sinks <- c(2,3,1,5)
socio <- electrograph(cbind(sources,sinks))

plot(socio)

}

