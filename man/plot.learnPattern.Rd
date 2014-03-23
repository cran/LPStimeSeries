\name{plot.learnPattern}
\alias{plot.learnPattern}
\title{Plot method for learnPattern objects}
\description{
  Plot the MSE of a learnPattern object over trees based on out-of-bag predictions
}
\usage{
\method{plot}{learnPattern}(x, type="l", main=deparse(substitute(x)), ...)
}
\arguments{
  \item{x}{an object of class \code{learnPattern}.}
  \item{type}{type of plot.}
  \item{main}{main title of the plot.}
  \item{...}{other graphical parameters.}
}
\value{
  Invisibly,  MSE of the \code{learnPattern} object.
}
\note{
  This function does not work for \code{learnPattern} if
  \code{oob.predict=FALSE} during training.

}
\seealso{
  \code{\link{learnPattern}}
}
\examples{
data(GunPoint)
ensemble=learnPattern(GunPoint$trainseries,oob.pred=TRUE,replace=TRUE)
plot(ensemble)
}
\author{Mustafa Gokce Baydogan}
\keyword{regression}
\keyword{tree}
