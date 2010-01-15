\name{pair.sequence}
\alias{pair.sequence}
\title{pair.sequence}
\description{
Given an input integer, produces a sequence of all n-choose-2 pairs in the data.
}
\usage{
pair.sequence(nn)
}
\arguments{
  \item{nn}{Maximum number of nodes.}
}
\value{
An (n-choose-2)-by-2 matrix whose rows take values (1,2), (1,3), ... (n-1,n).
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{
pair.sequence(5)
}
