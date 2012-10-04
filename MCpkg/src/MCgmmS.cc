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
#include "include/scythestat/optimize.h"
#include "MCfcds.h"
#include "MCrng.h"
#include "include/sandwich.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
#include <Rmath.h>
#undef pnorm // for name mangling with scythe: Rmath is C and does not know namespaces
#include <Rdefines.h>
#include <Rinternals.h>  // defines handling of R objects from C
#include <omp.h>



//using namespace std;


scythe::Matrix<> user_fun_gmms(SEXP fun, SEXP par, SEXP data, SEXP myframe) {

	SEXP R_fcall;
	if(!isFunction(fun)) error("`fun' must be a function");
	if(!isEnvironment(myframe)) error("`myframe' must be an environment");

	//PROTECT(R_fcall = lang3(fun, R_NilValue, R_NilValue));
    PROTECT(R_fcall = lang3(fun, par, data));
	//SETCADR(R_fcall, data);

	//SETCADR(R_fcall, par);

	SEXP funval;
	PROTECT(funval = eval(R_fcall, myframe));

	if (!isReal(funval)) error("`fun' must return double values");

	double* fv_data = REAL(funval); // should be nobs x nmom

	const unsigned int fv_nr = nrows(funval);

	const unsigned int fv_nc = ncols(funval);

	scythe::Matrix<> fv_m(fv_nr, fv_nc, fv_data);

	//if (fv_m == R_PosInf) error("`fun' returned +Inf");
	//if (R_IsNaN(fv_m) || R_IsNA(fv_m)) error("`fun' returned NaN or NA");
	UNPROTECT(2);
	return fv_m;

}

struct obj_fun_gmms {

	// save the relevant model objects inside of the functor
	scythe::Matrix<> moments_;
	scythe::Matrix<> moments_sum_;
	scythe::Matrix<> weights_;
	SEXP fun;
	SEXP par;
	SEXP data;
	SEXP env;

	double operator() (scythe::Matrix<> theta) {

		double fv;
		PROTECT(par = allocMatrix(REALSXP, theta.rows(), theta.cols()));
		for(unsigned int i = 0; i < theta.rows(); ++i) {
           REAL(par)[i] = theta(i, 0);
		}
		moments_ = user_fun_gmms(fun, par, data, env); // is (nobs x nmom)
        for(unsigned int i = 0; i < moments_.rows(); ++i) {
        	moments_sum_ += t(moments_(i, scythe::_));
        }
        //Rprintf("moments_sum_: %10.5f %10.5f %10.5f\n", moments_sum_(0,0), moments_sum_(1,0), moments_sum_(2,0));
        fv = (t(moments_sum_) * weights_ * moments_sum_)(0,0);
        //Rprintf("unscaled fv = %10.5f\n", fv);
        fv *= (1.0 / moments_.rows()) * 10e-5;
        //Rprintf("fv = %23.5f\n", fv);
        //Rprintf("Number of observations: %i\n", moments_.rows());
        /*Rprintf("Optimal par: \n\n");
        Rprintf("phi: %10.5f\n", theta(0,0));
        Rprintf("theta: %10.5f\n", theta(1,0));
        Rprintf("rho: %10.5f\n", theta(2,0));*/
        //Rprintf("-----------------------------\n");
        UNPROTECT(1);

        return fv;
	}
};

struct obj_fun_gmms2 {

	// save the relevant model objects inside the functor
	scythe::Matrix<> data_;
	scythe::Matrix<> moments_sum_;
	scythe::Matrix<> weights_;
	scythe::Matrix<> residuals_;
	double funv;


    double operator() (scythe::Matrix<> theta) {

        residuals_ = data_(scythe::_, 0) - (theta(0,0) + theta(1,0)) * data_(scythe::_, 1) + (theta(0,0) + theta(2,0) * theta(1,0)) * data_(scythe::_, 2);

        double moment1, moment2, moment3 = 0;
        #pragma omp parallel for reduction(+:moment1, moment2, moment3) schedule(dynamic)
        for(unsigned int i = 0; i < data_.rows(); ++i) {
        	moment1 += data_(i, 1) * data_(i, 2) - data_(i, 1) * data_(i, 1) * theta(2,0);
        	moment2 += residuals_(i) * data_(i, 1);
        	moment3 += residuals_(i) * data_(i, 2);
        }
        //Rprintf("moment1: %10.5f\t moment2: %10.5f\t moment3: %10.5f\n", moment1, moment2, moment3);
        moments_sum_(0,0) = moment1;
        moments_sum_(0,1) = moment2;
        moments_sum_(0,2) = moment3;
        funv = (moments_sum_ * weights_ * t(moments_sum_)) (0,0);

        funv *= 1.0 / data_.rows();
        return funv;

    }
};

struct moments_deriv {
	scythe::Matrix<> data_m;
	scythe::Matrix<> moment_deriv_m;

	scythe::Matrix<> operator() (const scythe::Matrix<>& theta_v, const unsigned int index) {

		  const unsigned int npar = theta_v.rows();
          scythe::Matrix<> moment1_d(1, npar);
          scythe::Matrix<> moment2_d(1, npar);
          scythe::Matrix<> moment3_d(1, npar);

          moment1_d(0,0) = 0.0;
          moment1_d(0,1) = 0.0;
          moment1_d(0,2) = -data_m(index, 1) * data_m(index, 1);

          moment2_d(0,0) = (-data_m(index, 1) + data_m(index, 2)) * data_m(index, 1);
          moment2_d(0,1) = (-data_m(index, 1) + theta_v(2,0) * data_m(index, 2)) * data_m(index, 1);
          moment2_d(0,2) = theta_v(1,0) * data_m(index, 2) * data_m(index, 1);

          moment3_d(0,0) = (-data_m(index, 1) + data_m(index, 2)) * data_m(index, 2);
          moment3_d(0,1) = (-data_m(index, 1) + theta_v(2,0) * data_m(index, 2)) * data_m(index, 2);
          moment3_d(0,2) = theta_v(1,0) * data_m(index, 2) * data_m(index, 2);

          moment_deriv_m = scythe::rbind(moment1_d, moment2_d);
          moment_deriv_m = scythe::rbind(moment_deriv_m, moment3_d);

          return moment_deriv_m;
	}
};

struct moments {

	scythe::Matrix<> data_m;
	scythe::Matrix<> moments_m;
    double tmp;
	scythe::Matrix<> operator() (const scythe::Matrix<>& theta_v) {

        #pragma omp parallel for schedule(dynamic)
	    for(unsigned int i = 0; i < data_m.rows(); ++i) {
		   moments_m(i, 0) = data_m(i, 1) * data_m(i, 2) - data_m(i, 1) * data_m(i, 1) * theta_v(2,0);
		   tmp = data_m(i, 0) - (theta_v(0,0) + theta_v(1,0)) * data_m(i, 1) - (theta_v(0,0) + theta_v(2,0) * theta_v(1,0)) * data_m(i, 2);
		   moments_m(i, 1) = tmp * data_m(i, 1);
		   moments_m(i, 2) = tmp * data_m(i, 2);

	    }

	    return moments_m;
	}
};

template<typename RNGTYPE>
void MCgmmS_impl(rng<RNGTYPE>& stream, SEXP& fun, SEXP& myframe,
		const unsigned int niter, const unsigned int nobs,
		const Matrix<>& par, const Matrix<>& covM, SEXP& sample_SEXP) {

	Matrix<> mu(3, 1);
	Matrix<> variables(nobs, 3);
	Matrix<> sample_copula(3,1);
	Matrix<> sample_u(3,1);

	for(unsigned int o = 0; o < nobs; ++o) {

        sample_copula = stream.rmvnorm(mu, covM);
        sample_u = scythe::pnorm(sample_copula, 0, 1);
        if(o == 0) {
        	variables(o, 0) = qbinom(sample_u(0,0), 1, 0.5, 1, 0);
        	if(variables(o, 0) == 0) variables(o, 0) = -1;
        }
        else {
        	if(sample_u(0,0) <= par(2,0)) variables(o, 0) = variables(o - 1, 0);
        	else variables(o, 0) = variables(o - 1, 0) * (-1);
        }
        variables(o, 1) = qnorm(sample_u(1,0), 0, 1, 1, 0);
        variables(o, 2) = qnorm(sample_u(2,0), 0, 1, 1, 0);
	}

	// construct the sample
	Matrix<> sample(100000, 4);
	sample(_, 1) = variables(0, 0, nobs - 2, 0);
	sample(_, 2) = variables(1, 0, nobs - 1, 0);
	sample(0, 3, nobs - 2, 3) = variables(0, 1, nobs - 2, 1) + variables(0, 2, nobs -2, 2) - variables(1, 2, nobs - 1, 2);
    sample(0, 0, nobs - 2, 0) = sample(0, 1, nobs - 2, 1) * (par(0,0) + par(1,0)) - sample(0, 2, nobs - 2, 2) * (par(0,0)
    		+ 0.2 * par(1,0)) + sample(0, 3, nobs - 2, 3);

    // put matrix into sample_SEXP
	for (unsigned int i = 0; i < (nobs - 1); ++i) {
		for (unsigned int j = 0; j < 4; ++j) {
		  REAL(sample_SEXP)[i + (nobs - 1) * j] = sample(i,j);
		}
	}

	// GMM two-step method:
	/* first step */
	scythe::Matrix<> init_par(par.rows(), par.cols());
	init_par = 0.05, 0.05, 0.1;
	scythe::Matrix<> weights_m = scythe::eye(par.rows());
	obj_fun_gmms2 obj_fun;
	obj_fun.residuals_ = scythe::Matrix<>(nobs, 1);
	obj_fun.data_ = sample;
	obj_fun.moments_sum_ = scythe::Matrix<>(1, par.rows());
	obj_fun.weights_ = weights_m;
	Matrix<> opt_par = scythe::BFGS(obj_fun, init_par, stream, 100, 10e-6, false);

	/* second step */
	moments moments;
	moments.data_m = sample;
	moments.moments_m = scythe::Matrix<>(sample.rows(), 3);
	scythe::Matrix<> moments_m = moments(opt_par);
	weights_m = sandwich::meat(moments_m, true);
	weights_m = scythe::invpd(weights_m);
	init_par = opt_par;
	obj_fun.weights_ = weights_m;
	opt_par = scythe::BFGS(obj_fun, init_par, stream, 100, 10e-6, false);

	/* results */
	Rprintf("Optimal par: \n\n");
	Rprintf("phi: %10.5f\n", opt_par(0,0));
	Rprintf("theta: %10.5f\n", opt_par(1,0));
	Rprintf("rho: %10.5f\n", opt_par(2,0));
	moments.moments_m = moments(opt_par);
	moments_deriv momderiv;
	momderiv.data_m = sample;
    scythe::Matrix<> neweyWestM = sandwich::NeweyWest(opt_par, moments_m, weights_m, momderiv);
    Rprintf("Newey West Cov:\n");
    Rprintf("%10.10f\t%10.10f\t%10.10f\n", neweyWestM(0, 0), neweyWestM(0, 1), neweyWestM(0,2));
    Rprintf("%10.10f\t%10.10f\t%10.10f\n", neweyWestM(1, 0), neweyWestM(1, 1), neweyWestM(1,2));
    Rprintf("%10.10f\t%10.10f\t%10.10f\n\n", neweyWestM(2, 0), neweyWestM(2, 1), neweyWestM(2,2));


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
 	   const unsigned int nobs_intern = nobs + 1;

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

       MCPKG_PASSRNG2MODEL(MCgmmS_impl, fun, myframe, niter, nobs_intern, par,
    		   covM, sample_SEXP);
       UNPROTECT(1);

       return sample_SEXP;
    }
}

#endif


