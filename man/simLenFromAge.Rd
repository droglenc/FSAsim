% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simLenFromAge.R
\name{simLenFromAge}
\alias{simLenFromAge}
\title{Simulate fish length from given ages.}
\usage{
simLenFromAge(ages, Linf, K, t0, sigma, additive = TRUE, digits = 0)
}
\arguments{
\item{ages}{A vector containg ages of individual fish.}

\item{Linf}{A numeric representing the asyptotic mean length (L_infinity) in
the von Bertalanffy growth model.}

\item{K}{A numeric representing the Brody growth coefficient (K) in the von
Bertalanffy growth model.}

\item{t0}{A numeric representing the time when the mean length is zero (t_0)
in the von Bertalanffy growth model.}

\item{sigma}{A numeric representing the standard deviation (i.e., individual
error) around the von Bertalanffy growth model.}

\item{additive}{A logical indicating whether the standard deviation is for
the additive- (\code{TRUE}; default) or multiplicative-error version of the
von Bertalanffy growth model.}

\item{digits}{A numeric controlling the number of digits to which the length
data should be rounded before returning.}
}
\value{
A vector containing the simulated lengths for individual fish.
}
\description{
Constructs simulated fish lengths from a set of given ages and parameters for
a von Bertalanffy growth model supplied by the user.
}
\details{
This simulation simulates fish lengths by first predicing the mean
length-at-age given a fish's age using the von Bertalanffy growth model and
then randomly selecting a length from a normal distribution given this mean
length-at-age and the supplied value for sigma.
}
\examples{
## Load FSA package for Summarize()
require(FSA)

## set seed for repeatability
set.seed(5234734)

## Simulated individual ages (random)
#    see simAges functions
bg.ages <- simAges(N0=500,A=0.35)

## Simulated lengths, given the above ages
bg.lens <- simLenFromAge(bg.ages,228,0.206,0,8)

## Some summaries
df <- data.frame(age=bg.ages,len=bg.lens)
Summarize(len~age,data=df,digits=1)
plot(len~age,data=df,pch=16,col=rgb(0,0,0,0.1),xlab="Age",ylab="Length")

}
\seealso{
\code{\link{simAges}}, \code{\link{simLenSelectP}},
\code{\link{simLenSelectM}}
}
\keyword{misc}
