\name{learnPattern}
\alias{learnPattern}
\alias{learnPattern.default}
\alias{print.learnPattern}

\title{Learn Local Auto-Patterns for Time Series Representation and Similarity}
\description{
  \code{learnPattern} implements ensemble of regression trees (based on
  Breiman and Cutler's original Fortran code) to learn local auto-patterns 
  for time series representation. Ensemble of regression trees are used to
  learn an autoregressive model. A local time-varying autoregressive 
  behavior is learned by the ensemble.
}
\usage{
\method{learnPattern}{default}(x,
   segment.factor=c(0.05,0.95),
   random.seg=TRUE, target.diff=TRUE, segment.diff=TRUE, 
   random.split=0,
   ntree=200,
   mtry=1,
   replace=FALSE,
   sampsize=if (replace) ceiling(0.632*nrow(x)) else nrow(x),
   maxdepth=6,
   nodesize=5,
   do.trace=FALSE,
   keep.forest=TRUE,
   oob.pred=FALSE,
   keep.errors=FALSE, 
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
  \item{random.split}{Type of the split. If set to zero (0), splits are generated based on decrease in SSE in target segment
	Setting of one (1) generates the split value randomly between max and min values. Setting of two (2) generates a kd-tree
	type of split (i.e. median of the values at each node is chosen as the split).}
  \item{ntree}{Number of trees to grow.  Larger number of trees are preferred
    	if there is no concern regarding the computation time.}
  \item{mtry}{Number of predictor segments randomly sampled as candidates at each
   	split.  Note that it is preset to 1 for now.}
  \item{replace}{Should bagging of time series be done with replacement? All training time series are used
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
  \item{oob.pred}{if replace is set to \code{TRUE}, predictions for the time
     series observations are returned. }   
  \item{keep.errors}{If set to \code{TRUE}, the mean square error (MSE) of target prediction over target segments is evaluated
   for each tree. If \code{oob.pred=TRUE}, this information is evaluated on ``out-of-bag'' samples at each tree. }
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
  \item{segment.length}{used segment length settings by the trees of ensemble}     
  \item{nobs}{number of observations in a segment}
  \item{ntree}{number of trees grown}
  \item{maxdepth}{maximum depth level for each tree}
  \item{mtry}{number of predictor segments sampled for spliting at each node.}
  \item{target}{starting time of the target segment for each tree.}
  \item{target.type}{type of the target segment; 1 if observed series, 2 if difference
	series.}  
  \item{forest}{a list that contains the entire forest; \code{NULL} if
    \code{keep.forest=FALSE}.}
  \item{oobprediction}{predicted observations based on ``out-of-bag'' time series 
    are returned if \code{oob.pred=TRUE}}
  \item{ooberrors}{Mean square error (MSE) over the trees evaluated using the predicted 
  observations on ``out-of-bag'' time series is returned if \code{oob.pred=TRUE}. }
  \item{inbag}{\code{n} by \code{ntree} matrix be returned that keeps track of 
    which samples are ``in-bag'' in which trees if \code{keep.inbag=TRUE}}
  \item{errors}{Mean square error (MSE) of target prediction over target segments
   for each tree. If \code{oob.pred=TRUE}, Mean square error (MSE) is reported based on ``out-of-bag'' samples at each tree.}
}
\note{
  OOB predictions may have missing values (i.e. \code{NA}) if time series is not
  left out-of-bag during computations. Even, it is left out-of-bag, there is a
  potential of some observations (i.e. time frames) not being selected as the target.
  In such cases, there will no OOB predictions.
}

\references{
  Baydogan, M. G. (2013), ``Learned Pattern Similarity``, Homepage: 
  \url{http://www.mustafabaydogan.com/learned-pattern-similarity-lps.html}.
  
  Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1),
  5-32.
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

## Learn patterns randomly on GunPoint training series with default parameters
ensemble=learnPattern(GunPoint$trainseries, random.split=1)

## Find the similarity between test and training series and classify test series
similarity=computeSimilarity(ensemble,GunPoint$testseries,GunPoint$trainseries)
NearestNeighbor=apply(similarity,1,which.min)
predicted=GunPoint$trainclass[NearestNeighbor]
accuracy=sum(predicted==GunPoint$testclass)/nrow(GunPoint$testseries)
print(100*accuracy)

## Learn patterns by training each tree on a random subsample
## and classify test time series
ensemble=learnPattern(GunPoint$trainseries,replace=TRUE)
similarity=computeSimilarity(ensemble,GunPoint$testseries,GunPoint$trainseries)
NearestNeighbor=apply(similarity,1,which.min)
predicted=GunPoint$trainclass[NearestNeighbor]
print(predicted)

## Learn patterns and do predictions on OOB time series
ensemble=learnPattern(GunPoint$trainseries,replace=TRUE,target.diff=FALSE,oob.pred=TRUE)
## Plot first series and its OOB approximation
plot(GunPoint$trainseries[1,],xlab='Time',ylab='Observation',
	type='l',lty=1,lwd=2)
points(c(1:ncol(GunPoint$trainseries)),ensemble$oobpredictions[1,],
	type='l',col=2,lty=2,lwd=2)
legend('topleft',c('Original series','Approximation'),
	col=c(1,2),lty=c(1,2),lwd=2)

}
\keyword{similarity}
\keyword{tree}
\keyword{classification}% at least one, from doc/KEYWORDS
\keyword{regression}

