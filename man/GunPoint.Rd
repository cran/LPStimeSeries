\name{GunPoint}
\alias{GunPoint}
\docType{data}
\title{The Gun-Point Data}
\description{
  This is the Gun-Point data from The UCR Time Series Database.
}
\usage{
data(GunPoint)
}
\format{
  \code{GunPoint} is a list with one training time series dataset and one test
  time series dataset provided as separate matrices. There are 50 cases (rows)
  for training dataset with 150 variables (columns). Similarly there are 150
  cases for test dataset with 150 variables. Variables are representing the
  observations over time. In other words, they are ordered so that a row is
  a univariate time series. Originally, this is a classification problem where
  there are two classes. Therefore, list stores the class information for both
  training and test time series. This information is stored in arrays of length 50
  and 150 for training and test time series respectively (so each time series
  is associated with a class).
  
  Description by Chotirat Ann Ratanamahatana and Eamonn Keogh in their publication
  ``Everything you know about Dynamic Time Warping is Wrong`` is as follows:
   
  ``...This dataset comes from the video surveillance domain. The dataset has two classes, 
  each containing 100 instances. All instances were created using one female actor and 
  one male actor in a single session. The two classes are: Gun-Draw: The actors have 
  their hands by their sides. They draw a replicate gun from a hip-mounted holster, 
  point it at a target for approximately one second, then return the gun to the holster, 
  and their hands to their sides. Point: The actors have their gun by their sides. 
  They point with their index fingers to a target for approximately one second, and 
  then return their hands to their sides. For both classes, we tracked the centroid 
  of the actor's right hands in both X- and Y-axes, which appear to be highly correlated; 
  therefore, in this experiment, we only consider the X-axis for simplicity...``  
}
\source{
   The original data is at \url{http://www.cs.ucr.edu/~eamonn/time_series_data/}.
}
\references{
  Ratanamahatana, C. A. and Keogh. E. (2004). Everything you know about Dynamic Time Warping is Wrong. 
  In proceedings of SIAM International Conference on Data Mining (SDM05),  pp.506-510
  Newport Beach, CA, April 21-23	
}
\seealso{
  \code{\link{learnPattern}}, \code{\link{computeSimilarity}}
}
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
print(predicted)
}
\author{Mustafa Gokce Baydogan}
\keyword{datasets}
