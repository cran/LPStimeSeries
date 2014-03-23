\name{computeSimilarity}
\alias{computeSimilarity}
\title{Compute similarity between time series based on learned patterns}
\description{
  Compute similarity between time series. Raw time series can be provided together
  with \code{learnPattern} object so that the representation for the time series are
  generated internally and similarity is computed based on these representations. The
  other option is to provide the representations (instead of raw time series) and to
  compute the similarity without a need for \code{learnPattern} object.
}
\usage{
computeSimilarity(object=NULL,testseries=NULL,refseries=NULL,
   maxdepth=NULL,which.tree=NULL,terminal=TRUE,testrepresentation,refrepresentation)
}
\arguments{
  \item{object}{an object of class \code{learnPattern}.}
  \item{refseries}{reference time series.}
  \item{testseries}{test time series.}
  \item{maxdepth}{maximum depth level to be used to generate representations
  for similarity computations.}
  \item{which.tree}{array of trees to be used for similarity computation.}
  \item{terminal}{\code{TRUE} if similarity is computed over the learned representations.}
  \item{testrepresentation}{learned representation for test time series.}
  \item{refrepresentation}{learned representation for reference time series.}
}
\value{
 A similarity matrix of size ``the number of test series`` by ``the number of reference
 series`` is returned. Similarity between test series and reference series is
 defined as the number of mismatching patterns based on the representation generated
 by the trees. See LPS paper for details. 
}

\note{
 Similarity matrix can also be computed over representations if it is generated 
 using \code{\link{predict.learnPattern}}. This will probably take longer time
 compared to computing the similarity directly using the ensemble. However, if you
 are using LPS for retrieval purposes, bounding schemes (such as early abondon) can
 be used (requires further implementation) with the learned representations.
}

\references{
  Baydogan, M. G. (2013), ``Learned Pattern Similarity``, Homepage: 
  \url{http://www.mustafabaydogan.com/learned-pattern-similarity-lps.html}.
}

\seealso{
  \code{\link{learnPattern}}, \code{\link{predict.learnPattern}}
}
\examples{

data(GunPoint)
set.seed(71)
## Learn patterns on GunPoint training series with default parameters
ensemble=learnPattern(GunPoint$trainseries)

## Find the similarity between test and training series
sim=computeSimilarity(ensemble,GunPoint$testseries,GunPoint$trainseries)

## Find similarity using representations, 
## First generate representations
trainRep=predict(ensemble, GunPoint$trainseries, nodes=TRUE)
testRep=predict(ensemble, GunPoint$testseries, nodes=TRUE)

## Then compute the similarity (city-block distance), 
## takes longer but we keep the representation
sim2=computeSimilarity(testrepresentation=testRep,refrepresentation=trainRep)

## Find the similarity based on first 100 trees
sim=computeSimilarity(ensemble,GunPoint$testseries,GunPoint$trainseries,which.tree=c(1:100))


}
\author{Mustafa Gokce Baydogan}
\keyword{similarity}
\keyword{tree}
