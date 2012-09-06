#ifndef MCGMM_CC
#define MCGMM_CC

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

#include <Rdefines.h>
#include <Rinternals.h>  // defines handling of R objects from C


using namespace std;
using namespace scythe;

double user_fun_gmm(SEXP fun, SEXP par, SEXP myframe) {
	
	SEXP R_fcall;
	if(!isFunction(fun)) error("`fun' must be a function");
	if(!isEnvironment(myframe)) error("myframe must be an environment");
	
	PROTECT(R_fcall = lang2(fun, R_NilValue));
	SETCADR(R_fcall, vartheta);
	SEXP funval;
	PROTECT(funval = eval(R_fcall, myframe));
	
	if (!isReal(funval)) error("`fun' must return a double");
	double fv = REAL(funval); //TODO: Make sure you put in the GMM objective function
	if (fv == R_PosInf) error("`fun' returned +Inf");
	if (R_IsNaN(fv) || R_IsNA(fv)) error("`fun' returned NaN or NA");
	UNPROTECT(2);
	return fv;
	
}

template<typename RNGTYPE>
void MCgmm_impl (rng<RNGTYPE>& stream, SEXP& fun,
		         SEXP& parlist, unsigned int nobs, unsigned int niter,
		         SEXP& myframe, bool verbose, const Matrix<>& covM,
		         SEXP& sample_SEXP)
{
	// define constants 
	const unsigned int nouter_par = length(parlist);
	const unsigned int npar = length(REAL(VECTOR_ELT(parlist,0))[0]);
	const Matrix<> covMc = cholesky(covM);

	// Initialize matrix to hold the sample
	Matrix<> sample(niter, npar, false);

	// Initialize list to hold all samples
	std::vector<Matrix<> > sampleL = std::vector(nouter_par);

	// put parameters into a scythe Matrix
	double* par_data = REAL(parlist);


}
extern "C" {

   SEXP MCgmm_cc(SEXP fun, SEXP myframe,
		   SEXP parameterList_R, SEXP arList_R, SEXP maList_R,
		   SEXP nobs_R, SEXP niter_R,
		   SEXP covM_R, SEXP verbose, SEXP ar_R, SEXP ma_R,
           SEXP dim_modCov_R,
           SEXP lecuyer_R, SEXP seedarray_R,
           SEXP lecuyerstream_R)
   {
	   // put rng stuff together
	   int seedarray[6];
	   for(int i=0; i<6; ++i) seedarray[i] = INTEGER(seedarray_R)[i];
	   int uselecuyer_cc = INTEGER(lecuyer_R)[0];
	   int lecuyerstream_cc = INTEGER(lecuyerstream_R)[0];
	   int* uselecuyer = &uselecuyer_cc;
	   int* lecuyerstream = &lecuyerstream_cc;
	   
	   // put covM_R into a Matrix
	   double* covM_data = REAL(covM_R);
	   const int covM_nr = nrows(covM_R);
	   const int covM_nc = ncols(covM_R);
	   Matrix <> covM (covM_nc, covM_nr, covM_data);
	   covM = t(covM);
	   
	   const unsigned int nmodels = length(parameterList_R);
	   const unsigned int niter = niter;
	   const unsigned int npar = length(VECTOR_ELT(parameterList_R,0));
	   
	   // set container to hold simulation values
	   SEXP sample_SEXP;
	   PROTECT(sample_SEXP = allocVector(VECSXP, nmodels)); //TODO: Consider to add diagnostics
	   for (unsigned int m = 0; m < nmodels; ++m) {
		   SEXP sample_vec_SEXP;
		   PROTECT(sample_vec_SEXP = allocMatrix(REALSXP, niter, npar));
           SET_VECTOR_ELT(sample_SEXP, m, sample_vec_SXP);
	   }
	   Rprintf("I am arrived");
	   
	   
   }
}
