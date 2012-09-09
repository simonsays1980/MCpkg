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
#include <Rmath.h>

#include <Rdefines.h>
#include <Rinternals.h>  // defines handling of R objects from C


using namespace std;
using namespace scythe;


double user_fun_gmm(SEXP fun, SEXP par, SEXP myframe) {
	
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

/*double user_fun_quantile(SEXP fun, SEXP args, SEXP myframe) {

	SEXP R_fcall;
	if(!isFunction(fun)) error("`fun' must be a function");
	if(!isEnvironment(myframe)) error("`myframe' must be an environment");

	SEXP docall;
	PROTECT(docall = mkChar("do.call"));
	PROTECT(R_fcall = lang3(docall, args, R_NilValue));
	SEXP funval;
	PROTECT(funval = eval(R_fcall, myframe));

	if(!isReal(funval)) error("`fun' must return a double");
	double fv = REAL(funval)[0];
	if (fv == R_PosInf || fv == R_NegInf) error("`fun' returned +/- Inf");
	if (R_IsNaN(fv) || R_IsNA(fv)) error("`fun' returned NaN or NA");
	UNPROTECT(3);
	return fv;

}*/

double user_fun_quantile(SEXP expr, SEXP val, SEXP args, SEXP env) {

	SEXP fnamec;
	unsigned int n = length(args);
	PROTECT(fnamec = expr = allocList(n + 2));
	SET_TYPEOF(fnamec, LANGSXP);
	SETCAR(fnamec, expr);
	fnamec = CDR(fnamec);
	SETCAR(fnamec, val);
	fnamec = CDR(fnamec);
	Rprintf("I am here: line 76 \n");
    for(unsigned int i = 0; i < n; ++i) {
    	SETCAR(fnamec, VECTOR_ELT(args, i));
    }
    Rprintf("I am here: line 80 \n");
    SEXP funval;
    SEXP call;
    PROTECT(funval = eval(fnamec, env));
    PROTECT(funval = eval(call, env));
    Rprintf("I am here: line 84 \n");

    double fv = REAL(call)[0];
    Rprintf("I am here: line 86 \n");
    if(fv == R_PosInf || fv == R_NegInf) error("'quantile function' returned +/- Inf");
    if(R_IsNaN(fv) || R_IsNA(fv)) error("quantile function' returned NaN or NA");
    UNPROTECT(3);
    return fv;

}

template<typename RNGTYPE>
void MCgmm_impl (rng<RNGTYPE>& stream, SEXP& fun, SEXP& myframe,
		         SEXP& parlist, SEXP& mreglist, SEXP& merrlist,
		         SEXP& mregparlist,
		         SEXP& arlist, SEXP& malist,
		         const unsigned int nobs, const unsigned int niter,
		         const unsigned int nar, const unsigned int nma,
		         const unsigned int verbose, const Matrix<>& covM,
		         SEXP& sample_SEXP)
{

	// define constants 
	const unsigned int nmodels = length(parlist);
	const unsigned int npar = length(VECTOR_ELT(parlist,0));
	const unsigned int nreg = length(VECTOR_ELT(mreglist, 0));
	const unsigned int nerr = length(VECTOR_ELT(merrlist, 0));
	const unsigned int ncopula_dim = nreg + nerr;

	// initialize matrix of means
	Matrix<> mu(ncopula_dim, 1);
	Rprintf("rows mu: %i \n", mu.rows());
	Rprintf("rows covM: %i \n", covM.rows());

	// THE MONTE CARLO SAMPLING
	for ( unsigned int m = 0; m < nmodels; ++m ) { // models

		// initialize matrix to hold the sample of estimated parameters
		Matrix<> sample_par(niter, npar, false);
		unsigned int narg = length(VECTOR_ELT(VECTOR_ELT(mregparlist, m), 0));

		// assign lists of margins and margin parameters
		SEXP sample_mrpList;
		PROTECT(sample_mrpList = allocVector(VECSXP, narg));
        sample_mrpList = VECTOR_ELT(mregparlist, m);
        SEXP sample_mrList;
        PROTECT(sample_mrList = allocVector(VECSXP, narg));
        sample_mrList = VECTOR_ELT(mreglist, m);

		for (unsigned int i = 0; i < niter; ++i) { // iterations

            // sample the copula
			Matrix<> sample_copV = stream.rmvnorm(mu, covM);

			Rprintf("rows sample_copV: %i \n", sample_copV.rows());
			Rprintf("cols sample_copV: %i \n", sample_copV.cols());
			Rprintf("row %i: %10.5f\n", 1, sample_copV(0,0));

            SEXP rval;
            PROTECT(rval = ScalarReal(sample_copV(0,0)));
            SEXP args;
            /*PROTECT(args = allocVector(REALSXP, narg));*/
            PROTECT(args = allocVector(VECSXP, narg));
            unsigned int nmpar = length(VECTOR_ELT(sample_mrpList, 0));
            for(unsigned int a = 0; a < nmpar; ++a) {
            	SET_VECTOR_ELT(args, a, VECTOR_ELT(sample_mrpList, 0))[a];
            	/*REAL(args)[a] = REAL(VECTOR_ELT(sample_mrpList, 0))[a];*/
            }
            Rprintf("I am here: line 148");
            SEXP qfun;
            PROTECT(qfun = VECTOR_ELT(mreglist, 0));
            double quantile = user_fun_quantile(qfun, rval, args, myframe);
            UNPROTECT(3);

		}
		UNPROTECT(2);
	}


}
extern "C" {

   SEXP MCgmm_cc(SEXP fun, SEXP myframe,
		   SEXP parameterList_R, SEXP margin_regList_R,
		   SEXP margin_reg_parList_R,
		   SEXP margin_errList_R, SEXP margin_err_parList_R,
		   SEXP copula, SEXP arList_R, SEXP maList_R,
		   SEXP nobs_R, SEXP niter_R,
		   SEXP covM_R, SEXP verbose, SEXP ar_R, SEXP ma_R,
           SEXP modeldim_R,
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
	   Matrix<> covM (covM_nc, covM_nr, covM_data);
	   
	   const unsigned int nmodels = length(parameterList_R);
	   const unsigned int niter = niter;
	   const unsigned int npar = length(VECTOR_ELT(parameterList_R,0));
	   
	   // set container to hold simulation values
	   SEXP sample_SEXP;
	   PROTECT(sample_SEXP = allocVector(VECSXP, nmodels)); //TODO: Consider to add diagnostics
	   for (unsigned int m = 0; m < nmodels; ++m) {
		   SEXP sample_vec_SEXP;
		   sample_vec_SEXP = allocMatrix(REALSXP, niter, npar); // TODO: Make sure all elements of the list are protected
           SET_VECTOR_ELT(sample_SEXP, m, sample_vec_SEXP);
	   }
	   
	   MCPKG_PASSRNG2MODEL(MCgmm_impl, fun, myframe, parameterList_R,
			   margin_regList_R, margin_errList_R, margin_reg_parList_R,
			   arList_R, maList_R, INTEGER(nobs_R)[0], INTEGER(niter_R)[0],
			   INTEGER(ar_R)[0], INTEGER(ma_R)[0], INTEGER(verbose)[0],
			   covM, sample_SEXP);
	   UNPROTECT(1);
       Rprintf("I am at the end");
	   // return the sample of parameters
       return sample_SEXP;

   }
}

#endif
