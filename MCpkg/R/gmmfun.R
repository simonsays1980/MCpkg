# TODO: Add comment
# 
# Author: simonzehnder
###############################################################################


"gmmfun" <- function(theta, data) {
	
	res = data[,1] - (theta[1] + theta[2]) * data[,2] + (theta[1] + theta[3] * theta[2]) * data[,3]
    
	moments = matrix(0, nrow=nrow(data), ncol=3)
	moments[,1] = data[,2] * data[,3] - data[,2]^2 * theta[3]
	moments[,2] = res * data[,2]
	moments[,3] = res * data[,3]

	return(moments)
}
