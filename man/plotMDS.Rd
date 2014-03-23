\name{plotMDS}
\alias{plotMDS}
\title{Multi-dimensional Scaling Plot of Learned Pattern Similarity}
\description{
  Plot the scaling coordinates of the Learned Pattern Similarity.
}
\usage{
plotMDS(object, newdata, classinfo, k=2, palette=NULL, pch=20, ...) 
}
\arguments{
  \item{object}{an object of class \code{learnPattern}, as that
    created by the function \code{learnPattern}.}
  \item{newdata}{a data frame or matrix containing the data for similarity
    computation.}
  \item{classinfo}{labels for the time series for color-coding.}
  \item{k}{number of dimensions for the scaling coordinates.}
  \item{palette}{colors to use to distinguish the classes; length must
    be the equal to the number of levels.}
  \item{pch}{plotting symbols to use.}
  \item{...}{other graphical parameters.}
}
\value{
  The output of \code{\link[stats]{cmdscale}} on scaled Learned Pattern
  similarity is returned invisibly.
}
\note{
  If \code{k > 2}, \code{\link{pairs}} is used to produce the
  scatterplot matrix of the coordinates.
  
  The entries of the similarity matrix is divided by the maximum possible
  similarity which is \code{2*sum(object$nobs)}
}
\seealso{\code{\link{learnPattern}}}
\examples{
set.seed(1)
data(GunPoint)
## Learn patterns on GunPoint training series with default parameters
ensemble=learnPattern(GunPoint$trainseries)
plotMDS(ensemble, GunPoint$trainseries,GunPoint$trainclass)

## Using different symbols for the classes:
plotMDS(ensemble, GunPoint$trainseries,GunPoint$trainclass, 
         palette=rep(1, 2), pch=as.numeric(GunPoint$trainclass))
         
## Learn patterns on GunPoint training series with random splits
ensemble=learnPattern(GunPoint$trainseries,random.split=1)
plotMDS(ensemble, GunPoint$trainseries,GunPoint$trainclass,main='Random Splits')

}
\author{Mustafa Gokce Baydogan}
\keyword{similarity}
