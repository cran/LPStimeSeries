\name{predict.learnPattern}
\alias{predict.learnPattern}
\title{predict method for \code{learnPattern} objects}
\description{
  Representation generation for test data using learnPattern.
}
\usage{
\method{predict}{learnPattern}(object, newdata, which.tree=FALSE,
   nodes=TRUE, maxdepth=NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{learnPattern}, as that
    created by the function \code{learnPattern}.}
  \item{newdata}{a data frame or matrix containing new data.}
  \item{which.tree}{\code{FALSE} if the representation is needed to be
  generated over all trees of ensemble. Set to an integer value if the representation
  is required to be generated for one tree specified by the value set.}
  \item{nodes}{\code{TRUE} generates the representation based on the trees.
  \code{TRUE} is the only option for now. A real-valued prediction is still under
  development. }
  \item{maxdepth}{The maximum depth level to generate the representation }
  \item{...}{not used currently.}
}

\value{
	Returns the learned pattern representation for the time series in the dataset.
	Basically, it is the count of observed patterns at each terminal node.
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
}
\keyword{similarity}
