# TODO: Add comment
# 
# Author: simonzehnder
###############################################################################


"MCgmmS" <- function(fun, parameters = NULL, niter = 1000, nobs = 1000,
		covM = NULL, seed = NA, ...) {
	
	# form seed
	seeds <- form.seeds(seed)
	lecuyer <- seeds[[1]]
	seed.array <- seeds[[2]]
	lecuyer.stream <- seeds[[3]]
	
	# setup the environment so that fun can see the things passed as ...
	gmm.fun <- function(ttt) fun(ttt, ...)
	env.gmm.fun <- environment(fun = gmm.fun)
	
	
	.Call("MCgmmS_cc", fun, env.gmm.fun, as.double(parameters), as.integer(niter), 
			as.integer(nobs), as.matrix(covM), as.integer(lecuyer), 
			as.integer(seed.array), as.integer(lecuyer.stream));
} 
