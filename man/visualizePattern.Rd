\name{visualizePattern}
\alias{visualizePattern}
\title{Plot of the patterns learned by the ensemble of the regression trees}
\description{
  \code{visualizePattern} visualizes the patterns implied by the terminal
  nodes of the trees from \code{learnPattern} object.
}
\usage{
visualizePattern(object, x, which.terminal, orient=c(2,2))
}
\arguments{
  \item{object}{an object of class \code{learnPattern}, as that
    created by the function \code{learnPattern}.}
  \item{x}{a data frame or matrix containing the data for pattern 
    visualization.}
  \item{which.terminal}{id of the terminal node determining the decision
    rules to be used for identifying patterns}
  \item{orient}{orientation of the plot (determines the grid structure and 
    how many patterns to be visualized).}
}
\value{
   A list with the following components are returned invisibly. 
  \item{predictor}{predictor segments residing in the \code{which.terminal}.}
  \item{target}{target segments implied by the \code{which.terminal}.}
  \item{tree}{the tree id corresponding to the \code{which.terminal}.}
  \item{terminal}{the id of the terminal node for the \code{tree}.}     
}
\note{
 Patterns are visualized for the time series for which the frequency of the 
 observations in the pattern is the largest. If more than one plot is requested
 through the setting of \code{orient}, the patterns are plotted for the time
 series based on the descending order of the frequency. 
 
 Currently, patterns are visualized based on the first predictor segment 
 (sampled at the root node). This visualization can be done based on the
 predictor segment sampled at each level of the tree.
 
 \code{predictor} and \code{target} are of size \code{x} where the patterns
 are numerical values and the rest of the entries are \code{NA}s. 
}
\seealso{\code{\link{learnPattern}},\code{\link{predict.learnPattern}}}
\examples{
set.seed(71)
data(GunPoint)
## Learn patterns on GunPoint training series with default parameters
ensemble=learnPattern(GunPoint$trainseries)

## Find representations
trainRep=predict(ensemble, GunPoint$trainseries, nodes=TRUE)

## Find the average frequency over the terminal nodes
avgFreq=apply(trainRep,2,mean)

## Find the terminal node that has the maximum average and visualize
termid=which.max(avgFreq)
visualizePattern(ensemble,GunPoint$trainseries,termid,c(2,1))

}
\author{Mustafa Gokce Baydogan}
\keyword{similarity}
