\name{getTreeInfo}
\alias{getTreeInfo}
\title{Extract a single tree from the ensemble.}
\description{
  This function extracts the structure of a tree from a
  \code{learnPattern} object.
}
\usage{
getTreeInfo(object, which.tree=1)
}
\arguments{
  \item{object}{a \code{\link{learnPattern}} object.}
  \item{which.tree}{which tree to extract?}
}
                                        
\value{
  is a list with the following components:
  \item{segment.length}{the proportion of the time series length used for 
     both predictors and targets.}
  \item{target}{starting time of the target segment.}
  \item{target.type}{type of the target segment; 1 if observed series, 2 
	if difference series.}
  \item{tree}{Tree structure matrix with seven columns and number of rows 
	equal to total number of nodes in the tree. }
	
  The seven columns of the \code{tree} structure matrix are:
  \item{left daughter}{the row where the left daughter node is; 0 if the
  node is terminal}
  \item{right daughter}{the row where the right daughter node is; 0 if
    the node is terminal}
  \item{split segment}{start time of the segment used to split the node}
  \item{split type}{type of the predictor segment used to split the node; 1 if 
  observed series, 2 if the different series are used. 0 if the node is 
  terminal}    
  \item{split point}{where the best split is}
  \item{status}{is the node terminal (-1) or not (-3)}
  \item{prediction}{the prediction for the node}
     
}
\note{
  For numerical predictors, data with values of the variable less than
  or equal to the splitting point go to the left daughter node.
}
%\references{
%}
\seealso{
  \code{\link{learnPattern}}
}
\examples{
data(GunPoint)
set.seed(71)

## Learn patterns on GunPoint training series with 50 trees
ensemble=learnPattern(GunPoint$trainseries,ntree=50)
getTreeInfo(ensemble, 3)
}
\author{Mustafa Gokce Baydogan}
\keyword{tree}
