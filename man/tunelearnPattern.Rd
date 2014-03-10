\name{tunelearnPattern}
\alias{tunelearnPattern}

\title{Tune Parameters of LPS for Time Series Classification}
\description{
  \code{tunelearnPattern} implements parameter selection for LPS in time series
  classification problems. LPS requires the setting of segment length (if segment
  length is not random) and depth parameter. Given training time series and alternative
  parameter settings, the best set of parameters that minimizes the cross-validation
  error rate is returned.
}
\usage{
tunelearnPattern(x, y, unlabeledx=NULL, nfolds=5, 
   segmentlevels=c(0.25,0.5,0.75), 
   mindepth=4, maxdepth=8, depthstep=2, 
   ntreeTry=25, target.diff=TRUE, segment.diff=TRUE, ...) 
}
\arguments{
  \item{x}{time series database as a matrix in UCR format. Rows are univariate time series, columns 
  are observations (for the \code{print} method, a \code{learnPattern} object).}
  \item{y}{labels for the time series given by \code{x}}
  \item{unlabeledx}{unlabeled time series dataset. Introduced for future purposes as LPS
  can benefit from unlabeled data.}
  \item{nfolds}{number of cross-validation folds for parameter evaluation.}
  \item{segmentlevels}{alternative segment level settings to be evaluated. Settings
  are provided as an array.}
  \item{mindepth}{minimum depth level to be evaluated.}
  \item{maxdepth}{maximum depth level to be evaluated.}
  \item{depthstep}{step size to determine the depth levels between \code{mindepth}
  and \code{maxdepth} to be evaluated.}
  \item{ntreeTry}{number of trees to be train for each fold.}
  \item{target.diff}{Can target segment be a difference feature?}
  \item{segment.diff}{Can predictor segments be difference feature?}  
  \item{...}{optional parameters to be passed to the low level function
    \code{tunelearnPattern}.}
}

\value{
  A list with the following components:
 
  \item{params}{evaluated parameter combinations as a matrix where rows are parameter 
  combinations and columns represent the settings. First and seconds columns are 
  the evaluated segment length level and depth respectively.}
  \item{errors}{cross-validation error rate  for each parameter combinations}
  \item{best.error}{the minimum cross-validation error rate obtained.}
  \item{best.seg}{the segment length level that provides the minimum cross-validation error.}
  \item{best.depth}{the depth level that provides the minimum cross-validation error.}
}

\references{
  Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1),
  5-32.

  Breiman, L (2002), ``Manual On Setting Up, Using, And Understanding
  Random Forests V3.1``, \url{http://oz.berkeley.edu/users/breiman/Using_random_forests_V3.1.pdf}.
  
  Baydogan, M. G. (2013), ``Learned Pattern Similarity``, Homepage: 
  \url{http://www.mustafabaydogan.com/learned-pattern-similarity-lps.html}.
}
\author{Mustafa Gokce Baydogan \email{baydoganmustafa@gmail.com}, 
 based on original Fortran code by Leo Breiman and Adele Cutler,
 R port by Andy Liaw and Matthew Wiener.}

\seealso{\code{\link{learnPattern}}, \code{\link{computeSimilarity}}}

\examples{
data(GunPoint)
set.seed(71)

##Tune segment length level and depth on GunPoint training series
tuned=tunelearnPattern(GunPoint$trainseries,GunPoint$trainclass)
print(tuned$best.error)
print(tuned$best.seg)
print(tuned$best.depth)

##Use tuned parameters to learn patterns
ensemble=learnPattern(GunPoint$trainseries,segment.factor=tuned$best.seg,
					  maxdepth=tuned$best.depth)

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
\keyword{tree}
\keyword{classification}
\keyword{regression}

