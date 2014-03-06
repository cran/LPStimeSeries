\name{learnPattern}
\alias{learnPattern}
\alias{learnPattern.default}
\alias{print.learnPattern}

\title{Learn Patterns for Time Series Representation and Similarity}
\description{
  \code{learnPattern} implements ensemble of regression trees (based on
  Breiman and Cutler's original Fortran code) to learn patterns for time series representation.
}
\usage{
\method{learnPattern}{default}(x,
   segment.factor=c(0.1,0.5),
   random.seg=TRUE, target.diff=TRUE, segment.diff=TRUE, 
   ntree=500,
   mtry=1,
   replace=FALSE,
   sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
   maxdepth=6,
   nodesize=5,
   do.trace=FALSE,
   keep.forest=TRUE,
   keep.inbag=FALSE, ...)
\method{print}{learnPattern}(x, ...)
}
\arguments{
  \item{x}{time series database as a matrix in UCR format. Rows are univariate time series, columns are
    observations (for the \code{print} method, a \code{learnPattern} object).}
  \item{segment.factor}{The proportion of the time series length to be used for 
     both predictors and targets, if \code{random.seg} is \code{TRUE} (default), minimum and maximum factor
     should be provided as array of length two.}
  \item{random.seg}{\code{TRUE} if segment length is random between thresholds 
	defined by \code{segment.factor}}
  \item{target.diff}{Can target segment be a difference feature?}
  \item{segment.diff}{Can predictor segments be difference feature?}
  \item{ntree}{Number of trees to grow.  Larger number of trees are preferred
    	if there is no concern regarding the computation time; }
  \item{mtry}{Number of predictor segments randomly sampled as candidates at each
   	split.  Note that it is preset to 1 for now.}
  \item{replace}{Should bagging of time series be done? All training time series are used
   	if \code{FALSE} (default).}
  \item{sampsize}{Size(s) of sample to draw with replacement if replace is set
	to \code{TRUE}}
  \item{maxdepth}{The maximum depth of the trees in the ensemble. }
  \item{nodesize}{Minimum size of terminal nodes.  Setting this number
    	larger causes smaller trees to be grown (and thus take less time).}
  \item{do.trace}{If set to \code{TRUE}, give a more verbose output as
    \code{learnPattern} is run.  If set to some integer, then running
    output is printed for every \code{do.trace} trees.}
  \item{keep.forest}{If set to \code{FALSE}, the forest will not be
    retained in the output object. }
  \item{keep.inbag}{Should an \code{n} by \code{ntree} matrix be
    returned that keeps track of which samples are ``in-bag'' in which
    trees}
  \item{...}{optional parameters to be passed to the low level function
    \code{learnPattern}.}
}

\value{
  An object of class \code{learnPattern}, which is a list with the
  following components:

  \item{call}{the original call to \code{learnPattern}.}
  \item{type}{\code{regression}}
  \item{segment.factor}{the proportion of the time series length to be used for 
     both predictors and targets.}
  \item{nobs}{number of observations in a segment}
  \item{ntree}{number of trees grown}
  \item{maxdepth}{maximum depth level for each tree}
  \item{mtry}{number of predictor segments sampled for spliting at each node.}
  \item{target}{starting time of the target segment for each tree.}
  \item{forest}{a list that contains the entire forest; \code{NULL} if
    \code{keep.forest=FALSE}.}
  \item{inbag}{\code{n} by \code{ntree} matrix be returned that keeps track of 
    which samples are ``in-bag'' in which trees if \code{keep.inbag=TRUE}}
}

\references{
  Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1),
  5-32.

  Breiman, L (2002), ``Manual On Setting Up, Using, And Understanding
  Random Forests V3.1``, \url{http://oz.berkeley.edu/users/breiman/Using_random_forests_V3.1.pdf}.
  
  Baydogan, M. G. (2013), ``Learned Pattern Similarity``,
  
  \url{http://www.mustafabaydogan.com/learned-pattern-similarity-lps.html}.
}
\author{Mustafa Gokce Baydogan \email{baydoganmustafa@gmail.com}, 
 based on original Fortran code by Leo Breiman and Adele Cutler,
 R port by Andy Liaw and Matthew Wiener.}

\seealso{\code{\link{predict.learnPattern}}, \code{\link{computeSimilarity}}, \code{\link{tunelearnPattern}}}

\examples{
data(GunPoint)
set.seed(71)

## Learn patterns on GunPoint training series with default parameters
ensemble=learnPattern(GunPoint$trainseries)
print(ensemble)

## Find the similarity between test and training series based on the learned model
similarity=computeSimilarity(ensemble,GunPoint$testseries,GunPoint$trainseries)

## Find the index of 1 nearest neighbor (1NN) training series for each test series
NearestNeighbor=apply(similarity,1,which.min)

## Predicted class for each test series
predicted=GunPoint$trainclass[NearestNeighbor]

## Compute the percentage of accurate predictions
accuracy=sum(predicted==GunPoint$testclass)/nrow(GunPoint$testseries)
print(100*accuracy)

}
\keyword{similarity}
\keyword{tree}
\keyword{classification}% at least one, from doc/KEYWORDS
\keyword{regression}

