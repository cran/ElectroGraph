\name{electro.frats}
\docType{data}
\alias{electro.frats}
\alias{bernard.killworth.b}
\alias{newcomb}
\title{Two Fraternity Data Sets}
\description{
Contains the time series from the Newcomb fraternity study, as well as the 
observational set from the Bernard-Killworth human interaction data sets.
}
\usage{electro.frats}
\format{
newcomb: an array, 17-by-17-by-15. Each of the 15 array slices represents a time point; each row represents the social preferences of each individual for the others at the time. Originally rank choices of 1-16, these are now fractions (17 - (1-16))/17 such that higher preferences have higher connection.

bernard.killworth.b: a 1934-by-3 array. Columns 1 and 2 ("P1","P2") represent the dyad under observation; column 3 ("count") is the number of times these individuals were seen to interact over the observational period.
}
\source{
Originally obtained from the UCINET IV data source; included in R format for convenience and demonstration.

http://vlado.fmf.uni-lj.si/pub/networks/data/Ucinet/UciData.htm
}
\references{
Bernard H, Killworth P and Sailer L. (1980). Informant accuracy in social network data IV. Social Networks, 2, 191-218.
Newcomb T. (1961). The acquaintance process. New York: Holt, Reinhard & Winston.
}
\keyword{datasets}
