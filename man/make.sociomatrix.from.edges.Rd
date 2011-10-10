\name{make.sociomatrix.from.edges}
\alias{make.sociomatrix.from.edges}
\title{make.sociomatrix.from.edges}
\description{
Given an edge list, with or without values, produce the corresponding sociomatrix.
}
\usage{
make.sociomatrix.from.edges(inputmat, size=NULL, symmetric=FALSE, fidelities=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputmat}{A matrix with dimension n-by-(2,3,4); see below for
    details.}
  \item{size}{If non-null, specifies the number of nodes that should be
    in the resulting sociomatrix and implies that existing node labels
    are integers from 1 to size.}
  \item{symmetric}{If TRUE, creates a symmetric sociomatrix with undirected edges.}
  \item{fidelities}{If set, allows for ties that are perceived to be
    less than ``friends''. (Use with caution -- results have not been peer reviewed!}
}
\details{
  The first two columns of the matrix are the source and target of each directed edge (undirected if symmetry flag is set to TRUE).
  If inputmat has 3 columns, the third column represents the value of each tie.
  If inputmat has 4 columns, the third and fourth column represent the edge weights in each direction, as in (source->target), (target->source).
}
\value{
  A sociomatrix, where the number of rows and columns is equal to the count of unique nodes specified; row and column names are equal to the node names given.
}
\author{Andrew C. Thomas <act@acthomas.ca>}

\examples{
src <- 1:10
tgt <- 2:11
socio <- make.sociomatrix.from.edges(cbind(src,tgt), symmetric=TRUE)

values <- rgamma(10,3,3)
socio.value <- make.sociomatrix.from.edges(cbind(src,tgt,values), symmetric=TRUE)

back.values <- rgamma(10,3,12)
socio.value <- make.sociomatrix.from.edges(cbind(src,tgt,values,back.values))

}

