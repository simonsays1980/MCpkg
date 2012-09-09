# return name of the calling function
"calling.function" <-
  function(parentheses=TRUE) {
    calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
    if (parentheses){
      calling.function <- paste(calling.function, "()", sep="")
    }
    return(calling.function)
  }

# parse the passed seeds
# 1] if a scalar is passed, it is used by Mersennse twister
# 2] if a list of length two is passed, a parallel-friendly stream is
#    created using L'Ecuyer
"form.seeds" <-
  function(seed) {
    if(length(seed)==1) {
      if(is.na(seed)) seed <- 12345
      seed <- as.integer(seed)
      if(seed < 0) {
        cat("Error: Mersenne seed negative.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)                       
      }
      seeds <- list(0, rep(seed,6), 0)
    }
    if(length(seed)==2) {
      if(!is.list(seed)) {
        cat("Error: List must be passed to use L'Ecuyer.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      lec.seed <- seed[[1]]
      lec.substream <- as.integer(seed[[2]])
      if(is.na(lec.seed[1])) lec.seed <- rep(12345, 6)
      if(length(lec.seed) != 6) {
        cat("Error: L'Ecuyer seed not of length six.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if(!all(lec.seed >= 0))  {
        cat("Error: At least one L'Ecuyer seed negative.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if( max(lec.seed[1:3]) >= 4294967087){
        cat("Error: At least one of first three L'Ecuyer seeds\n")
        cat("  greater than or equal to 4294967087\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if( all(lec.seed[1:3] == 0 )){
        cat("Error: first three L'Ecuyer seeds == 0\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if( max(lec.seed[4:6]) >= 4294944443){
        cat("Error: At least one of last three L'Ecuyer seeds\n")
        cat("  greater than or equal to 4294944443\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }         
      if( all(lec.seed[4:6] == 0 )){
        cat("Error: last three L'Ecuyer seeds == 0\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if(lec.substream < 1) {
        cat("Error: L'Ecuyer substream number not positive.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)               
      }
      seeds <- list(1, lec.seed, lec.substream) 
    }
    if(length(seed)>2) {
      cat("Error: Seed passed as length greater than two.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)        
    }
    return(seeds)
  }

"check.mc.parameters" <- function(nobs, niter) {
  if( !is.numeric(nobs) ) {
    cat("Error: MC observations must be of type numeric. \n")
    stop("Respecify and call ", calling.function(), " agan. \n", 
         call. = FALSE)
  }
  if( !is.numeric(niter) ) {
    cat("Error: MC iterations not numeric. \n")
    stop("Please respecify and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  if( nobs < 0) {
    cat("Error: MC observations negative. \n")
    stop("Please respecify and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  if( niter < 0) {
    cat("Error: MC iterations negative. \n")
    stop("Please respecify and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  
}

"check.mc.input" <- function(fun = NULL, formula = NULL, ar.list = NULL, ma.list = NULL, 
                             parameter.list = NULL, margin.regressor.list = NULL, 
							 margin.regressor.parameter.list = NULL, copula = NULL,
							 margin.error.list = NULL, margin.error.parameter.list = NULL) {
						 
  sum.parameter <- ifelse( !is.null(parameter.list), sum(sapply(parameter.list, length)), 0)
  sum.ar <- ifelse( !is.null(ar.list), sum(sapply(ar.list, length)), 0)
  sum.ma <- ifelse( !is.null(ma.list), sum(sapply(ma.list, length)), 0)
  sum.mreg <- ifelse( !is.null(margin.regressor.list), sum(sapply(margin.regressor.list, length)), 0)
  sum.merr <- ifelse( !is.null(margin.error.list), sum(sapply(margin.error.list, length)), 0)
  sum.mreg.par <- ifelse( !is.null(margin.error.parameter.list), sum(sapply(margin.error.parameter.list, length)), 0)
  sum.merr.par <- ifelse( !is.null(margin.regressor.parameter.list), sum(sapply(margin.regressor.parameter.list, length)), 0)
  
  model.dim <- sum.ar + sum.ma + sum.mreg + sum.merr
  
  if ( is.null(fun) && is.null(formula) ) {
    cat("Error: no function specified. \n")
    stop("Respecify function and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  if ( is.null(parameter.list) ) {
    cat("Error: no or to less parameters set. \n")
    stop("Respecify parameters and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  
  if( !all(sapply(parameter.list, is.numeric)) ) {
    cat("Error: parameters are not of type numeric. \n")
    stop("Check input and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  
  # parameter vector lengths must be equal to the number of regressors
  if( sum.parameter != model.dim) {
    cat("Error: number of regressors and number of parameters differ. \n")
    warning("Check model and call ", calling.function(), " again. \n")
  } 
  
  # check if ar and ma lists contain only whole numbers
  r.ar <- sapply(ar.list, round)
  if( !is.null(ar.list) && !isTRUE(all.equal(ar.list, r.ar)) ) {
    cat("Error: autoregressive process wrongly specified. \n")
    stop("ar.list must be a list of whole numbers. \nCheck ar.list and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  r.ma <- sapply(ma.list, round)
  if( !is.null(ma.list) && !isTRUE(all.equal(ma.list, r.ma)) ) {
    cat("Error: moving average process wrongly specified. \n")
    stop("ma.list must be a list of whole numbers. \nCheck ma.list and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  
  if ( is.null(margin.error.list) ) {
	  cat("Error: at least one generating function for an error term has to be specified. \n")
	  stop("Specify a random number function and call ", calling.function(), " again\n", 
			  call. = FALSE)
  }
  
  # check if margin lists are correctly specified
  choices <- c("qnorm", "qunif", "qgamma", "qbeta", "qlnorm", "qchisq", 
		  "qnchisq", "qf", "qt", "qbinom", "qcauchy", "qexp", "qgeom",
		  "qhyper", "qnbinom", "qpois", "qweibull", "qlogis", "qnbeta",
		  "qnf", "qnt", "qtukey", "qwilcox", "qsignrank")
  res <- NULL
  try(res <- mapply(function(i) sapply(margin.error.list[[i]], match.arg, choices), 
				  seq(along=margin.error.list)), silent = TRUE)
  
  if( is.null(res) ) {
	  cat("Error: generating function does not exist. \n")
	  stop("at least one marginal function of margin.error.list does not exist. \nCheck random number generators and call ", calling.function(), " again. \n",
			  .call = FALSE)
  }
  
  if( is.null(margin.error.list) ) {
	  cat("Error: at least one error term has to be specified. \n")
	  stop("specify error term(s) and call ", calling.function(), " again. \n", 
			  call. = FALSE)
  }
  
  res <- NULL
  try(res <- mapply(function(i) sapply(margin.regressor.list[[i]], match.arg, choices), 
				  seq(along=margin.regressor.list)), silent = TRUE)
  
  if( is.null(res) ) {
    cat("Error: generating function does not exist. \n")
    stop("at least one marginal function of margin.regressor.list, does not exist. \nCheck random number generators and call ", calling.function(), " again. \n",
         call. = FALSE)
  } 
  if ( is.null(margin.regressor.list) ) {
	  cat("Warning: the simulation will be run with error term(s) only. \n")
	  warnings("this could result in unpredictable behavior of user-specified functions.\n")
  }
  
  
}