"learnPattern.default" <-
    function(x,
    	     segment.factor=c(0.1,0.5),
	         random.seg=TRUE, target.diff=TRUE, segment.diff=TRUE, 
             ntree=500,
             mtry=1,
             replace=FALSE,
             sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
             maxdepth = 6,
             nodesize = 5,
	         do.trace=FALSE,
             keep.forest=TRUE,
             keep.inbag=FALSE, ...) {
              
    n <- nrow(x)
    p <- ncol(x)

	if(length(segment.factor)>1){
		random.seg <- TRUE
		segment.factor <- sort(segment.factor)
	} else {
		random.seg <- FALSE
	}
    if (n == 0) stop("data (x) has 0 rows")
    
	ncat <- rep(1, p) #indicator for categorical variables (future purposes)
	
    ## overcome R's lazy evaluation:
    keep.forest <- keep.forest

    ## Make sure mtry is in reasonable range.
    if (mtry != 1)
        warning("invalid mtry: reset to within valid range")
    mtry <- 1

    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")

    ## Compiled code expects variables in rows and time series in columns.
    x <- t(x)
    storage.mode(x) <- "double"
    nt <- if (keep.forest) ntree else 1
    
    ## possible total # of nodes
    nrnodes <- 2^(maxdepth+1) - 1
	
	if(!replace){
		keep.inbag <- FALSE
	} else {
		sampsize <- ceiling(.632*n)
	}
		
	if(random.seg){
		segment.length <- runif(ntree,segment.factor[1],segment.factor[2])
	} else {
		segment.length <- rep(segment.factor[1],ntree)
	}

	rfout <- .C("regRF_time_series",
                    x,
                    as.double(segment.length),
                    as.integer(target.diff),
					as.integer(segment.diff),
                    as.integer(c(n, p)),
                    as.integer(sampsize),
                    as.integer(nodesize),
                    as.integer(nrnodes),
                    as.integer(ntree),
                    as.integer(mtry),
                    as.integer(ncat),
                    as.integer(do.trace),  
                    target = integer(ntree),           
                    targetType = integer(ntree),   
                    ndbigtree = integer(ntree),
                    nodedepth = matrix(integer(nrnodes * nt), ncol=nt),
                    nodestatus = matrix(integer(nrnodes * nt), ncol=nt),
                    splitType = matrix(integer(nrnodes * nt), ncol=nt),
                    leftDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                    rightDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                    nodepred = matrix(double(nrnodes * nt), ncol=nt),
                    bestvar = matrix(integer(nrnodes * nt), ncol=nt),
                    xbestsplit = matrix(double(nrnodes * nt), ncol=nt),
                    keep = as.integer(c(keep.forest, keep.inbag)),
                    replace = as.integer(replace),
                    inbag = if (keep.inbag)
                    matrix(integer(n * ntree), n) else integer(1),
                    PACKAGE="LPStimeSeries")[c(13:26)]
        ## Format the forest component, if present.
        if (keep.forest) {
            max.nodes <- max(rfout$ndbigtree)
            rfout$nodestatus <-
                rfout$nodestatus[1:max.nodes, , drop=FALSE]
            rfout$nodedepth <-
                rfout$nodedepth[1:max.nodes, , drop=FALSE]                
            rfout$splitType <-
                rfout$splitType[1:max.nodes, , drop=FALSE]               
            rfout$bestvar <-
                rfout$bestvar[1:max.nodes, , drop=FALSE]
            rfout$nodepred <-
                rfout$nodepred[1:max.nodes, , drop=FALSE]
            rfout$xbestsplit <-
                rfout$xbestsplit[1:max.nodes, , drop=FALSE]
            rfout$leftDaughter <-
                rfout$leftDaughter[1:max.nodes, , drop=FALSE]
            rfout$rightDaughter <-
                rfout$rightDaughter[1:max.nodes, , drop=FALSE]
        }
        cl <- match.call()
        cl[[1]] <- as.name("learnPattern")

        out <- list(call = cl,
                    type = "regression",
					segment.factor = segment.factor,
					segment.length = segment.length,
					nobs = floor(segment.length*p),
                    ntree = ntree,
                    maxdepth = maxdepth,
                    mtry = mtry,
                    target = rfout$target,
                    targetType = rfout$targetType,
                    forest = if (keep.forest)
                    c(rfout[c("ndbigtree", "nodedepth", "nodestatus", "splitType", "leftDaughter",
                              "rightDaughter", "nodepred", "bestvar",
                              "xbestsplit")],
                    list(nrnodes=max.nodes)) else NULL,
                    inbag = if (keep.inbag)
                    matrix(rfout$inbag, nrow(rfout$inbag),ntree) else NULL)
                           
    class(out) <- "learnPattern"
    return(out)
}
