\name{predict.learnPattern}
\alias{predict.learnPattern}
\title{predict method for \code{learnPattern} objects}
\description{
  Representation generation for test data using learnPattern.
}
\usage{
\method{predict}{learnPattern}(object, newdata, which.tree=NULL,
   nodes=TRUE, maxdepth=NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{learnPattern}, as that
    created by the function \code{learnPattern}.}
  \item{newdata}{a data frame or matrix containing new data.}
  \item{which.tree}{\code{NULL} if the representation is needed to be
  generated over all trees of ensemble. Set to an integer value if the representation
  is required to be generated for one tree specified by the value set.}
  \item{nodes}{\code{TRUE} generates the representation based on the trees.
.  \code{FALSE} generates a real-valued prediction for each time point.  }
  \item{maxdepth}{The maximum depth level to generate the representation }
  \item{...}{not used currently.}
}

\value{
	Returns the learned pattern representation for the time series in the dataset
	if \code{nodes} is set \code{TRUE}. Basically, it is the count of observed patterns at
	each terminal node. Otherwise predicted values for each time series in \code{newdata}
	are returned.
}
\references{
  Baydogan, M. G. (2013), ``Learned Pattern Similarity``, Homepage: 
  \url{http://www.mustafabaydogan.com/learned-pattern-similarity-lps.html}.
  
  Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1),
  5-32.
}
\author{Mustafa Gokce Baydogan}

\seealso{\code{\link{learnPattern}}}

\examples{
data(GunPoint)
set.seed(71)
## Learn patterns on GunPoint training series with default parameters
ensemble=learnPattern(GunPoint$trainseries)

## Find representations
trainRep=predict(ensemble, GunPoint$trainseries, nodes=TRUE)
testRep=predict(ensemble, GunPoint$testseries, nodes=TRUE)

## Check size of the representation for training data
print(dim(trainRep))

## Learn patterns on GunPoint training series (target cannot be difference series)
ensemble=learnPattern(GunPoint$trainseries,target.diff=FALSE)

## Predict observations for test time series
predicted=predict(ensemble,GunPoint$testseries,nodes=FALSE)

## Plot an example test time series 
plot(GunPoint$testseries[5,],type='l',lty=1,xlab='Time',ylab='Observation',lwd=2)
points(c(1:ncol(GunPoint$testseries)),predicted$predictions[5,],type='l',col=2,lty=2,lwd=2)
legend('topleft',c('Original series','Approximation'),col=c(1,2),lty=c(1,2),lwd=2)

}
\keyword{similarity}
