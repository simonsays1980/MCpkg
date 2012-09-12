/*
 * MCgmmS.cc
 *
 *  Created on: Sep 9, 2012
 *      Author: simonzehnder
 */

#ifndef MCGMMS_CC
#define MCGMMS_CC

#include "include/scythestat/matrix.h"
#include "include/scythestat/distributions.h"
#include "include/scythestat/stat.h"
#include "include/scythestat/la.h"
#include "include/scythestat/ide.h"
#include "include/scythestat/smath.h"
#include "MCfcds.h"
#include "MCrng.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
#include <Rmath.h>
#undef pnorm
#include <Rdefines.h>
#include <Rinternals.h>  // defines handling of R objects from C


//using namespace std;


double user_fun_gmms(SEXP fun, SEXP par, SEXP myframe) {

	SEXP R_fcall;
	if(!isFunction(fun)) error("`fun' must be a function");
	if(!isEnvironment(myframe)) error("`myframe' must be an environment");

	PROTECT(R_fcall = lang2(fun, R_NilValue));
	SETCADR(R_fcall, par);
	SEXP funval;
	PROTECT(funval = eval(R_fcall, myframe));

	if (!isReal(funval)) error("`fun' must return a double");
	double fv = REAL(funval)[0]; // TODO: Make sure you put in the GMM objective function
	if (fv == R_PosInf) error("`fun' returned +Inf");
	if (R_IsNaN(fv) || R_IsNA(fv)) error("`fun' returned NaN or NA");
	UNPROTECT(2);
	return fv;

}

template<typename RNGTYPE>
void MCgmmS_impl(rng<RNGTYPE>& stream, SEXP& fun, SEXP& myframe,
		const unsigned int niter, const unsigned int nobs,
		const Matrix<>& par, const Matrix<>& covM, SEXP& sample_SEXP) {

	Matrix<> mu(2, 1);
	Matrix<> variables(nobs, 4);
	Matrix<> sample_copula(2,1);
	Matrix<> sample_u(2,1);

    #pragma omp parallel for
	for(int o = 0; o < nobs; ++o) {

        sample_copula = stream.rmvnorm(mu, covM);
        sample_u = scythe::pnorm(sample_copula, 0, 1);
        variables(o, 1) = qbinom(sample_u(0,0), 1, 0.5, 1, 0);
        if(variables(o, 1) == 0) variables(o, 1) = -1;
        variables(o, 2) = qnorm(sample_u(1,0), 0, 1, 1, 0);
	}

	// put matrix into sample_SEXP
	for (unsigned int i = 0; i < nobs; ++i) {
	    for (unsigned int j = 0; j < 4; ++j) {
	      REAL(sample_SEXP)[i + nobs * j] = variables(i,j);
	    }
	}

}
extern "C" {

   SEXP MCgmmS_cc(SEXP fun, SEXP myframe, SEXP parameters_R, SEXP niter_R,
    		SEXP nobs_R, SEXP covM_R,
    		SEXP lecuyer_R, SEXP seedarray_R, SEXP lecuyerstream_R) {

    	// put rng stuff together
 	   int seedarray[6];
 	   for(int i=0; i<6; ++i) seedarray[i] = INTEGER(seedarray_R)[i];
 	   int uselecuyer_cc = INTEGER(lecuyer_R)[0];
 	   int lecuyerstream_cc = INTEGER(lecuyerstream_R)[0];
 	   int* uselecuyer = &uselecuyer_cc;
 	   int* lecuyerstream = &lecuyerstream_cc;

 	   const unsigned int npar = length(parameters_R);
 	   const unsigned int niter = INTEGER(niter_R)[0];
 	   const unsigned int nobs = INTEGER(nobs_R)[0];

 	    // put covM_R into a Matrix
       double* covM_data = REAL(covM_R);
       const int covM_nr = nrows(covM_R);
       const int covM_nc = ncols(covM_R);
       Matrix<> covM (covM_nr, covM_nc, covM_data);

       // put parameters_R into a Matrix
       double* par_data = REAL(parameters_R);
       const int par_nr = length(parameters_R);
       const int par_nc = 1;
       Matrix<> par (par_nr, par_nc, par_data);

       SEXP sample_SEXP;
       PROTECT(sample_SEXP = allocMatrix(REALSXP, nobs, 4));

       MCPKG_PASSRNG2MODEL(MCgmmS_impl, fun, myframe, niter, nobs, par,
    		   covM, sample_SEXP);
       UNPROTECT(1);

       return sample_SEXP;
    }
}

#endif


