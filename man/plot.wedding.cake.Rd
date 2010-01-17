\name{plot.wedding.cake}
\alias{plot.wedding.cake}
\title{plot.wedding.cake}
\description{
Given an electrograph object for a valued relational data set, create a series of plots that illustrate the graph at various threshold levels for tie existence.
}
\usage{
plot.wedding.cake(x, connectivity.mode="sociomatrix", force.mode="fruchterman.reingold",
                  ego.focus=NULL, filebase="electrograph", #type="png",
                  lower.bound=NULL, upper.bound=NULL, main.title=TRUE, 
		  plot.width = 600, plot.height = 600, ...) 


}
\arguments{
  \item{x}{An object of class ``electrograph''.}
  \item{connectivity.mode}{Should the projected distances be calculated using the shortest path, or the equivalent conductance of the electro-social approach?}
  \item{force.mode}{Which force-directed placement method should be used to find node coordinates?}
  \item{ego.focus}{A vector of nodes whose importance should be considered highest in the plot.}
  \item{filebase}{An identifier for the series of images produced in the approach.}
  \item{lower.bound, upper.bound}{The minimum and maximum values for an edge to be plotted.}
  \item{main.title}{If true, print the bound values in the title of the plot.}
  \item{plot.width, plot.height}{The height and width in pixels of the output PNG files.}
  \item{...}{Additional options to pass to ``plot''.}
}
\value{
	Creates a series of plots in the PNG format corresponding to the choices of bounds in a directory named for ``filebase''.	
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{

#A sample Erdos-Renyi type graph for wedding caking.
pts <- pair.sequence(10)
edges <- exp(rnorm(dim(pts)[1],0,1))
socio <- electrograph(cbind(pts,edges))

plot.wedding.cake(socio)

}

