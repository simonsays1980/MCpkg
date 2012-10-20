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
#undef pnorm // for name mangling with scythe: Rmath is C and does not know namespaces but only definitions
#include <Rdefines.h>
#include <Rinternals.h>  // defines handling of R objects from C
#include <omp.h>

scythe::Matrix<> user_fun_gmmsR(SEXP fun, SEXP par, SEXP data, SEXP myframe) {

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

struct obj_fun_gmmsR {

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
		moments_ = user_fun_gmmsR(fun, par, data, env); // is (nobs x nmom)
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
	double funv;
	double moment1, moment2, moment3;


    double operator() (scythe::Matrix<> theta_v) {
        double moment1 = 0; 
        double moment2 = 0; 
        double moment3 = 0;
        double tmp;
    	unsigned int nobs = data_.rows();
        omp_set_num_threads(2);
        #pragma omp parallel for reduction(+:moment1, moment2, moment3) private(tmp) schedule(dynamic)
        for(unsigned int i = 0; i < nobs; ++i) {
        	tmp = 0;
        	moment1 += data_(i, 2) * data_(i, 1) - data_(i, 1) * data_(i, 1) * theta_v(2,0);
        	tmp = data_(i, 0) - (theta_v(0,0) + theta_v(1,0)) * data_(i, 1) + (theta_v(0,0) + theta_v(2,0) * theta_v(1,0)) * data_(i, 2);
        	moment2 += tmp * data_(i, 1);
                moment3 += tmp * data_(i, 2);
        }
        //Rprintf("moment1: %10.5f\t moment2: %10.5f\t moment3: %10.5f\n", moment1, moment2, moment3);
        moments_sum_(0,0) = moment1;
        moments_sum_(0,1) = moment2;
        moments_sum_(0,2) = moment3;
        funv = (moments_sum_ * weights_ * t(moments_sum_)) (0,0);

        funv *= 1.0 / nobs;
        return funv;

    }
};

struct moments_deriv {
	scythe::Matrix<> data_m;

	scythe::Matrix<> operator() (const scythe::Matrix<>& theta_v, const unsigned int index) {
        scythe::Matrix<> moment_deriv_m(3,3);
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
	scythe::Matrix<> operator() (const scythe::Matrix<>& theta_v, const unsigned int nmom) {
        
        const unsigned int nobs = data_m.rows();
        double tmp = 0;
        scythe::Matrix<> moments_m(nobs, nmom);
        
        #pragma omp parallel for schedule(dynamic) shared(moments_m) private(tmp)
	    for(unsigned int i = 0; i < nobs; ++i) {
	       tmp = 0;
		   moments_m(i, 0) = data_m(i, 2) * data_m(i, 1) - data_m(i, 1) * data_m(i, 1) * theta_v(2,0);
		   tmp = data_m(i, 0) - (theta_v(0,0) + theta_v(1,0)) * data_m(i, 1) + (theta_v(0,0) + theta_v(2,0) * theta_v(1,0)) * data_m(i, 2);
		   moments_m(i, 1) = tmp * data_m(i, 1);
		   moments_m(i, 2) = tmp * data_m(i, 2);

	    }

	    return moments_m;
	}
};

template<typename RNGTYPE>
void MCgmmS_impl(scythe::rng<RNGTYPE>& stream, SEXP& fun, SEXP& myframe,
		const unsigned int niter, const unsigned int nobs_intern,
		scythe::Matrix<>& par, const scythe::Matrix<>& covM, unsigned int verbose, SEXP& sample_SEXP, SEXP& parsample_SEXP) {

    /* define constants */
	const unsigned int nobs = nobs_intern - 1;

	/* gnerate sample container */
	scythe::Matrix<> pmatrix(niter, 16);

	/* determine number of processors for OpenMP loop */
	if(niter <= 50) {
            Rprintf("\nOpenMP will use 20 processors\n\n");
            omp_set_num_threads(20);
        }
	else { 
            int nP = omp_get_num_procs();
	    Rprintf("\nOpenMP will use %i processors\n\n", nP);
            omp_set_num_threads(nP);
        }

	/* compute Markov chain transition probability */
    const double ptrans = (1 - par(2,0) * 1)/2;

    double dveps, dvres, dSSE, dSST, dwald;

    /* generate pseudo-random number generator */
    mersenne mc_rng;
    mc_rng.initialize(1);

    /*lecuyer mc_rng;
    unsigned long start_seed_array[6] = {1,2,3,4,5,6};
    mc_rng.SetPackageSeed(start_seed_array);*/

    /* Monte Carlo sampler */
    /*omp_set_num_threads(nP);
    #pragma omp parallel for schedule(dynamic) shared(pmatrix) private(mc_rng, dveps, dres, dSSE, dSST, dwald)
*/	for(unsigned int iter = 0; iter < niter; ++iter) {

		/* set matrices */
		scythe::Matrix<> mu(3, 1);
		scythe::Matrix<> variables(nobs_intern, 3);
		scythe::Matrix<> sample_copula(3,1);
		scythe::Matrix<> sample_u(3,1);
		scythe::Matrix<> dmid_m(nobs, 1);
		scythe::Matrix<> explained_m(nobs, 1);
		scythe::Matrix<> mean_m(1,3);

		/* reset the pseudo-random number generator for no overlap */
		/*mc_rng.initialize(1 + iter * nobs * 3);*/
		/*if(iter != 0) {
        unsigned long rseed = 12345 + iter * nobs * 3;
        unsigned long rseed_array [6];
        for(unsigned int i = 0; i < 6; ++i) {
        	rseed_array[i] = rseed;
        }
        mc_rng.SetSeed(rseed_array);
		}*/

		for(unsigned int o = 0; o < nobs_intern; ++o) {

			sample_copula = mc_rng.rmvnorm(mu, covM);
			sample_u = scythe::pnorm(sample_copula, 0, 1);
			if(o == 0) {
				variables(o, 0) = qbinom(sample_u(0,0), 1, 0.5, 1, 0);
				if(variables(o, 0) == 0) variables(o, 0) = -1;
			}
			else {
				if(sample_u(0,0) <= (1 - ptrans)) variables(o, 0) = variables(o - 1, 0);
				else variables(o, 0) = variables(o - 1, 0) * (-1);
			}
			variables(o, 1) = qnorm(sample_u(1,0), 0, par(3,0), 1, 0);
			variables(o, 2) = qnorm(sample_u(2,0), 0, par(4,0), 1, 0);
		}

		/* construct the sample */
		Matrix<> sample(nobs, 4);
		sample(_, 1) = variables(0, 0, nobs_intern - 2, 0);
		sample(_, 2) = variables(1, 0, nobs_intern - 1, 0);
		// TODO: CHeck if it works also with scythe::_
		sample(0, 3, nobs_intern - 2, 3) = variables(0, 1, nobs_intern - 2, 1) + variables(0, 2, nobs_intern -2, 2) - variables(1, 2, nobs_intern - 1, 2);
		sample(0, 0, nobs_intern - 2, 0) = sample(0, 1, nobs_intern - 2, 1) * (par(0,0) + par(1,0)) - sample(0, 2, nobs_intern - 2, 2) * (par(0,0)
				+ par(2,0) * par(1,0)) + sample(0, 3, nobs_intern - 2, 3);

		/* construct price and midquote difference */
        dmid_m = par(1,0) * sample(scythe::_, 2) + variables(0, 1, nobs - 1, 1);

		/* put matrix into sample_SEXP */
		/*if(iter == 0) {
			for (unsigned int i = 0; i < nobs; ++i) {
				for (unsigned int j = 0; j < 4; ++j) {
				  REAL(sample_SEXP)[i + nobs * j] = sample(i,j);
				}
			}
		}*/

		// GMM two-step method:
		/* first step */
		scythe::Matrix<> init_par(3, 1);
		init_par = 0.05, 0.05, 0.1;
		scythe::Matrix<> weights_m = scythe::eye(3);
		obj_fun_gmms2 obj_fun;
		obj_fun.data_ = sample;
		obj_fun.moments_sum_ = scythe::Matrix<>(1, 3);
		obj_fun.weights_ = weights_m;
		Matrix<> opt_par = scythe::BFGS(obj_fun, init_par, mc_rng, 100, 10e-6, false);

		/* second step */
		moments moments;
		moments.data_m = sample;
		scythe::Matrix<> moments_m = moments(opt_par, 3);
		weights_m = sandwich::meat(moments_m, true);
		weights_m = scythe::invpd(weights_m);
		init_par = opt_par.copy();
		obj_fun.weights_ = weights_m;
		opt_par = scythe::BFGS(obj_fun, init_par, mc_rng, 100, 10e-6, false);

		dveps = scythe::var(dmid_m - opt_par(1,0) * sample(scythe::_,2));
		explained_m = ((opt_par(0,0) + opt_par(1,0)) * sample(scythe::_, 1) - (opt_par(0,0) + opt_par(2,0) * opt_par(1,0)) * sample(scythe::_, 2));
		dvres = scythe::var(sample(scythe::_, 0) - explained_m);
		dSSE = nobs * scythe::var(explained_m);
		dSST = nobs * scythe::var(sample(scythe::_, 0));


			pmatrix(iter, 0) = nobs;
			pmatrix(iter, 1, iter, 5) = par(scythe::_, 0);
			pmatrix(iter, 6, iter, 8) = opt_par(scythe::_,0);
			pmatrix(iter, 12) = dvres;
			pmatrix(iter, 13) = dveps;
			pmatrix(iter, 14) = dSSE/dSST;


		/* results */
		if (iter % verbose == 0) {
		  Rprintf("Monte Carlo iteration %i of %i \n", (iter+1), niter);
		  Rprintf("parameter = \n");

			  for (unsigned int i = 0; i < 3; ++i)
			  Rprintf("%10.5f\n", pmatrix(iter, i + 6));

		}
		moments_m = moments(opt_par, 3);
		moments_deriv momderiv;
		momderiv.data_m = sample;
		scythe::Matrix<> neweyWestM = sandwich::NeweyWest(opt_par, moments_m, weights_m, momderiv);
			pmatrix(iter, 9, iter, 11) = scythe::sqrt(scythe::diag(neweyWestM))(scythe::_, 0);


		/*Rprintf("Newey West Cov:\n\n");
		Rprintf("%10.10f\t%10.10f\t%10.10f\n", neweyWestM(0, 0), neweyWestM(0, 1), neweyWestM(0,2));
		Rprintf("%10.10f\t%10.10f\t%10.10f\n", neweyWestM(1, 0), neweyWestM(1, 1), neweyWestM(1,2));
		Rprintf("%10.10f\t%10.10f\t%10.10f\n\n", neweyWestM(2, 0), neweyWestM(2, 1), neweyWestM(2,2));*/

		dwald = (t(opt_par) * invpd(neweyWestM) * opt_par)(0,0);
		dwald *= nobs;

			pmatrix(iter, 15) = dwald;

		R_CheckUserInterrupt();

	}

	/* SEXP-matrix is always column major ordered */
	for(unsigned int i = 0; i < niter; ++i) {
		for(unsigned int j = 0; j < 16; ++j) {
			REAL(parsample_SEXP)[i + niter * j] = pmatrix(i, j);
		}
	}

	Rprintf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");
	Rprintf("\tMonte Carlo simulation:\n"
			"\titerations:%i\t\t observations:%i\n\n", niter, nobs);
	Rprintf("\tparameter means are: \n");
	Rprintf("\tphi:\t%-10.5f\n"
			"\ttheta:\t%-10.5f\n"
			"\trho:\t%-10.5f\n\n", scythe::mean(pmatrix(scythe::_, 6)),
			scythe::mean(pmatrix(scythe::_, 7)), scythe::mean(pmatrix(scythe::_, 8)));
	Rprintf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");



}
extern "C" {

   SEXP MCgmmS_cc(SEXP fun, SEXP myframe, SEXP parameters_R, SEXP niter_R,
    		SEXP nobs_R, SEXP covM_R, SEXP verbose,
    		SEXP lecuyer_R, SEXP seedarray_R, SEXP lecuyerstream_R) {

    	/* put rng stuff together */
 	   int seedarray[6];
 	   for(int i=0; i<6; ++i) seedarray[i] = INTEGER(seedarray_R)[i];
 	   int uselecuyer_cc = INTEGER(lecuyer_R)[0];
 	   int lecuyerstream_cc = INTEGER(lecuyerstream_R)[0];
 	   int* uselecuyer = &uselecuyer_cc;
 	   int* lecuyerstream = &lecuyerstream_cc;

 	   /* set constant parameters */
 	   const unsigned int niter = INTEGER(niter_R)[0];
 	   const unsigned int nobs = INTEGER(nobs_R)[0];
 	   const unsigned int nobs_intern = nobs + 1;

 	    /* put covM_R into a scythe::Matrix */
       double* covM_data = REAL(covM_R);
       const int covM_nr = nrows(covM_R);
       const int covM_nc = ncols(covM_R);
       Matrix<> covM (covM_nr, covM_nc, covM_data);

       /* put parameters_R into a scythe::Matrix */
       double* par_data = REAL(parameters_R);
       const int par_nr = length(parameters_R);
       const int par_nc = 1;
       Matrix<> par (par_nr, par_nc, par_data);

       SEXP sample_SEXP;
       PROTECT(sample_SEXP = allocMatrix(REALSXP, nobs, 4));

       SEXP parsample_SEXP;
       PROTECT(parsample_SEXP = allocMatrix(REALSXP, niter, 16));

       MCPKG_PASSRNG2MODEL(MCgmmS_impl, fun, myframe, niter, nobs_intern, par,
    		   covM, INTEGER(verbose)[0], sample_SEXP, parsample_SEXP);

       UNPROTECT(2);

       return parsample_SEXP;
    }
}

#endif


