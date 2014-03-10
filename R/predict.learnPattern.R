"predict.learnPattern" <-
    function (object, newdata, which.tree=FALSE,
				nodes=TRUE, maxdepth=NULL, ...)
{
    if (!inherits(object, "learnPattern"))
        stop("object not of class learnPattern")
    if (is.null(object$forest)) stop("No forest component in the object")

    x <- newdata
    if (nrow(x) == 0)
        stop("newdata has 0 rows")
    if (any(is.na(x)))
        stop("missing values in newdata")
	
    if(is.null(maxdepth)) maxdepth <- object$maxdepth
    
    if(maxdepth>object$maxdepth) {
		maxdepth <- object$maxdepth
		warning("invalid depth: reset to the maximum depth provided during training!")
    }
		
    keep <- 1:nrow(x)
    rn <- rownames(x)
    if (is.null(rn)) rn <- keep

    mdim <- ncol(x)
    ntest <- nrow(x)
    nclass <- object$forest$nclass
    nrnodes <- object$forest$nrnodes
    
    ## get rid of warning:
    op <- options(warn=-1)
    on.exit(options(op))
    x <- t(data.matrix(x))

	if (nodes){
		keepIndex <- c("nodeRep","lenRep")
		if(which.tree>0){
			nodexts <- integer(ntest * object$forest$nrnodes)
		} else {
			nodexts <- integer(ntest * object$forest$nrnodes * object$ntree )
		}			
		ans <- .C("regForest_represent",  
				as.double(x),
				as.integer(ntest),
				as.integer(which.tree),
				as.double(object$segment.length),
				as.integer(mdim),
				as.integer(object$ntree),
				object$forest$leftDaughter,
				object$forest$rightDaughter,
				object$forest$nodestatus,
				object$forest$nodedepth,
				object$forest$nrnodes,
				object$forest$xbestsplit,
				object$forest$bestvar,
				object$forest$splitType,
				object$forest$ndbigtree,
				as.integer(maxdepth),
				nodeRep = nodexts,
				lenRep = integer(1),
				PACKAGE = "LPStimeSeries")[keepIndex]
				
		res=t(matrix(ans$nodeRep[1:(ans$lenRep*ntest)], nrow=ans$lenRep))
		
	} 
    res
} 
    
				   		                    
