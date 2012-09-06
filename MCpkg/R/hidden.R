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

"check.mc.input" <- function(fun = NULL, formula = NULL, ar.list, ma.list, 
                             parameter.list, nregressors, margin.error.list, 
                             margin.regressor.list) {
  
  
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
  if( is.null(nregressors) ) {
    cat("Error: no number of regressors specified. \n")
    warning("Check model and call ", calling.function(), " again. \n")
  } 
  if( !is.null(ar.list) && all(sapply(lapply(ar.list, is.element, c(0,1,-1)), all)) ) {
    cat("Error: Autoregressive process wrongly specified. \n")
    stop("ar must be a vector of zeros and/or +/- 1. \nCheck ar and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  if( !is.null(ma.list) && all(sapply(lapply(ma.list, is.element, c(0,1,-1)), all)) ) {
    cat("Error: Moving average process wrongly specified. \n")
    stop("ma must be a vector of zeros and/or +/- 1. \nCheck ma and call ", calling.function(), " again. \n",
         call. = FALSE)
  }
  if ( is.null(margin.error.list) ) {
	  cat("Error: at least one generating function for an error term has to be specified. \n")
	  stop("Specify a random number function and call ", calling.function(), " agqain\n", 
			  call. = FALSE)
  }
  
  choices <- c("qnorm", "qunif", "qgamma", "qbeta", "qlnorm", "qchisq", 
		  "qnchisq", "qf", "qt", "qbinom", "qcauchy", "qexp", "qgeom",
		  "qhyper", "qnbinom", "qpois", "qweibull", "qlogis", "qnbeta",
		  "qnf", "qnt", "qtukey", "qwilcox", "qsignrank")
  res <- NULL
  try(res <- sapply(margin.error.list, match.arg, choices), silent = TRUE)
  
  if( is.null(res) ) {
    cat("Error: Generating function does not exist. \n")
    stop("At least one marginal function of margin.error.list, does not exist. \nCheck random number generators and call ", calling.function(), " again. \n",
         call. = FALSE)
  } 
  if ( is.null(margin.regressor.list) ) {
	  cat("Warning: the simulation will be run with error term(s) only. \n")
	  warnings("GMM function could behave strangly.\n")
  }
  
  res <- NULL
  try(res <- sapply(margin.regressor.list, match.arg, choices), silent = TRUE)
  
  if( is.null(res) ) {
    cat("Error: Generating function does not exist. \n")
    stop("At least one marginal function of margin.regressor.list does not exist. \nCheck random number generators and call ", calling.function(), " again. \n",
         .call = FALSE)
  }
}