"MCgmm" <- function(fun = NULL, nregressors = NULL, 
		            nerrors = 1, 
                    ar.list = NULL, ma.list = NULL,
                    parameter.list = NULL, 
                    margin.regressor.list = NULL, 
                    copula = c("rmvnorm", "rmvt"),
                    margin.error.list = NULL,
                    model.covM = NULL,
                    verbose = 1,
                    force.samp = FALSE,
					seed = NA,
                    nobs = 1000, niter = 1000, ...) {
  
  # check mc input 
  check.mc.input(fun = fun, formula = NULL, ar.list = ar.list, ma.list = ma.list, parameter.list = parameter.list, 
		  nregressor = nregressors, margin.error.list = margin.error.list, margin.regressor.list = margin.regressor.list)
  
  # check mc parameters
  check.mc.parameters(nobs = nobs, niter = niter)
  
  # ma and ar terms
  n.ar <- ifelse( !is.null(ar.list), sum(sapply(ar.list, sum)), 0 )
  n.ma <- ifelse( !is.null(ma.list), sum(sapply(ma.list, sum)), 0 )
  dim.model.covM <- nregressors - n.ar + nerrors - n.ma
  
  # form seed
  seeds <- form.seeds(seed)
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
  # setup the environment so that fun can see the things passed as ...
  gmm.fun <- function(ttt) fun(ttt, ...)
  env.gmm.fun <- environment(fun = gmm.fun)
  
  # covariance matrix
  if( is.null(model.covM) ) {
    model.covM <- diag(1, dim.model.covM)
  }
  else { # model.covM provided by user
    if( nrow(model.covM) != ncol(model.covM) || nrow(model.covM) != dim.model.covM) {
      stop("model.covM is not correctly specified. \nSpecify model.covM and call MCgmm() again. \n",
           call. = FALSE)
    }
    if( !isSymmetric(model.covM) ) {
      stop("model.covM must be symmetric.  \nSpecify model.covM and call MCgmm() again. \n",
           call. = FALSE)
    }
    if( !all(is.element(diag(model.covM), c(1))) ) {
      warning("model.covM must be a correlation matrix. \nMatrix has been normalized. \n")
      d <- diag(model.covM)
      D <- diag(1/sqrt(d))
      model.covM <- D %*% model.covM %*% D
      if(verbose == TRUE) {
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
        cat("           Normalized model.covM: \n")
        print(model.covM)
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
      }
    }
    if( max(model.covM) > max(diag(model.covM)) ) {
      stop("Correlations in model.covM are greater than one. \nRespecify correlations and call MCgmm() again. \n",
           call. = FALSE)
    }
    
    # get choleski
    CC <- NULL
    try(CC <- chol(model.covM), silent = TRUE)
    cov.flag = 0
    if( is.null(CC) ) {
      
      if( force.samp == TRUE) {
        ## enforce positive definiteness
        require(Matrix)
        model.covM.new <- nearPD(model.covM)
        
        try(CC <- chol(model.covM.new), silent = TRUE)
        model.covM <- model.covM.new
        cov.flag = 1
      }
      else { # force.samp = FALSE
        cov.flag = 2
      }
      
      if( cov.flag == 1 ) {
        warning("model.covM not positive definite.\nSampling proceeded after enforcing
              positive definiteness. \n")
        if( verbose == TRUE) {
          cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
          cat("           Nearest positive definite model.covM: \n")
          print(model.covM.new)
          cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
        }
      }
      if( cov.flag == 2) {
        cat("covM not positive definite.\n")
        cat("Sampling (as.specified) cannot proceed.\n")
        stop("Check model.covM and call MClm() again. \n", 
             call. = FALSE)
      }
    }
  }
  
  # call C++ function to do the sampling 
  sample <- .Call("MCgmm_cc", gmm.fun, env.gmm.fun, 
                  as.list(parameter.list), as.list(margin.regressor.list),
				  as.list(margin.error.list),
				  as.list(ar.list), as.list(ma.list),
                  as.integer(nobs), as.integer(niter),
                  as.matrix(model.covM), as.integer(verbose), 
                  as.integer(n.ar), as.integer(n.ma),
                  as.integer(dim.model.covM),
                  lecuyer = as.integer(lecuyer), 
                  seed.array = as.integer(seed.array),
                  lecuyer.stream = as.integer(lecuyer.stream))
}