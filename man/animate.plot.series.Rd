\name{animate.plot.series}
\alias{animate.plot.series}
\title{animate.plot.series}
\description{
Create a series of plots that animate the transitions between subsequent states of a relational data system.
}
\usage{
animate.plot.series(plots.to.animate, intermediates=10, filebase="electro-animate",
                                plot.width=600, plot.height=600, ...)
}
\arguments{
  \item{plots.to.animate}{A list of objects of class "electrograph.plot", the output from a plot(electrograph) command.}
  \item{intermediates}{The number of frames to intersperse between source plots.}
  \item{filebase}{The name of the directory into which the output files will be saved.}
  \item{plot.width, plot.height}{The height and width in pixels of the output PNG files.}
  \item{...}{Additional arguments to pass to ``plot''.}
}
\value{
	Creates a series of plots in the PNG format in a directory named for ``filebase''.	Also produces a shell script to produce an animated GIF of the sequence using the program ImageMagick.
	Function outputs the list of files produced in the animation.
}
\author{Andrew C. Thomas <act@acthomas.ca>}
\examples{

#uses data "newcomb" in this set.
data(electro.frats)
newcomb.electro <- list(NA)
for (kk in 1:dim(newcomb)[3]) {
  newcomb.electro[[kk]] <- electrograph(newcomb[,,kk])
}

newcomb.plot <- list(NA)
newcomb.plot[[1]] <- plot(newcomb.electro[[1]],distance.mode="electro.social", just=TRUE)
for (kk in 2:length(newcomb.electro))
  newcomb.plot[[kk]] <- plot(newcomb.electro[[kk]],distance.mode="electro.social", 
										just=TRUE, previous.elec=newcomb.plot[[kk-1]])

animate.plot.series(newcomb.plot)

}

