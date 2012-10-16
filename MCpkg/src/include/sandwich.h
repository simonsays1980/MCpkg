/*
 * sandwich.h
 *
 *  Created on: Sep 16, 2012
 *      Author: simonzehnder
 */

#ifndef SANDWICH_H_
#define SANDWICH_H_

#include <scythestat/matrix.h>
#include <scythestat/defs.h>
#include <scythestat/ide.h>
#include <scythestat/la.h>
#include <scythestat/optimize.h>
#include <cmath>
#include <omp.h>
#include <R.h>
#include <Rdefines.h>

namespace sandwich {
    namespace {
       typedef unsigned int uint;
    }

    /* forward declarations */
    scythe::Matrix<> computeWeights(const double&, const uint&);
    template<typename FUNCTOR>
    scythe::Matrix<> bread(const scythe::Matrix<>&, const scythe::Matrix<>&, const scythe::Matrix<>&, FUNCTOR);
    scythe::Matrix<> meat(const scythe::Matrix<>&, const bool);

    template<typename FUNCTOR>
    inline scythe::Matrix<>
    NeweyWest(const scythe::Matrix<>& theta_v, const scythe::Matrix<>& moments_m, const scythe::Matrix<>& weights_m, FUNCTOR momderiv_f) {

    	const uint nobs = moments_m.rows();
    	const uint npar = theta_v.rows();
    	scythe::Matrix<> neweyWest(npar, npar);
        scythe::Matrix<> bread_m = bread(theta_v, moments_m, weights_m, momderiv_f);
    	neweyWest =  invpd(t(bread_m) * scythe::invpd(meat(moments_m, true)) * bread_m);


    	return neweyWest * 1.0/nobs;
    }

    template<typename FUNCTOR>
    inline scythe::Matrix<>
    bread(const scythe::Matrix<>& theta_v, const scythe::Matrix<>& moments_m, const scythe::Matrix<>& weights_m, FUNCTOR momderiv_f) {

    	const uint nobs = moments_m.rows();
    	const uint nmom = moments_m.cols();
    	const uint npar = theta_v.rows();
    	scythe::Matrix<> bread_m(nmom, npar);
    	scythe::Matrix<> momderiv_m(nmom, npar);
    	scythe::Matrix<> tmp_m(nmom, npar);

        /* Fastest implementation so far. Set number of threads to the number of available cores.
         * TODO: check what happens if the loop gets nested in the parallel loop of the main function!
         */
    	omp_set_num_threads(2);
        #pragma omp parallel
    	{
             #pragma omp for schedule(dynamic) firstprivate(momderiv_f)
    		 for(unsigned int i = 0; i < nobs; ++i) {
    			 tmp_m += momderiv_f(theta_v, i);

    		 }
             #pragma omp single nowait
			 {
				 momderiv_m += tmp_m;
			 }
    	}
		/*for(unsigned int j = 0; j < nobs; ++j) {
			momderiv_m += momderiv_f(theta_v, j);
		}*/

		momderiv_m*= 1.0/(nobs - 3);
		// This is the bread for GMM with weights matrix. In the optimal case it reduces to the one below.
		// bread_m = weights_m * momderiv_m * scythe::invpd(t(momderiv_m) * weights_m * momderiv_m);
        bread_m = momderiv_m;

        return bread_m;
    }

    inline scythe::Matrix<>
    meat(const scythe::Matrix<>& moments_m, const bool centered) {

    	const uint nobs = moments_m.rows();
    	const uint nmom = moments_m.cols();

    	double bandwidth = std::pow((double)nobs, 0.2);
        double stepsize = 1.0 / (1.0 + bandwidth);
    	uint pos_w = 1 + std::floor(1.0 / stepsize);

    	scythe::Matrix<> weights = computeWeights(stepsize, pos_w);

    	scythe::Matrix<> meat_m(nmom, nmom);
    	scythe::Matrix<> demean_m(moments_m);
    	if(centered) {
    		scythe::Matrix<> means = scythe::meanc(moments_m);
    		scythe::Matrix<> ones = scythe::ones(nobs, 1);
    		demean_m = moments_m - ones * means;
    	}
    	meat_m = t(demean_m) * demean_m * 0.5 * weights(0,0);
    	scythe::Matrix<> zero_m(nmom, nmom);

    	if(weights.rows() > 1) {
    		for(uint i = 1; i < weights.cols(); ++i) {
    			meat_m += weights(i, 0) * t(demean_m(0, 0, nobs - i - 1, nmom - 1)) * demean_m(i, 0, nobs - 1, nmom - 1);
    		}
    	}

    	meat_m += t(meat_m);
    	meat_m *= 1.0 / nobs;

    	return meat_m;
    }

    /* Implementation of Bartlett Kernel for Newey-West. Other HACs or HCs use other kernels not implemented yet. */
    inline scythe::Matrix<>
    computeWeights(const double& stepsize, const uint& pos_w) {

        scythe::Matrix<> weights(pos_w, 1);
        double tmp;
        for(unsigned int i = 0; i < pos_w; ++i) {
            tmp = 1.0 - i * pos_w;
            weights(i,0) = tmp;
        }

        return weights;
    }

}


#endif /* SANDWICH_H_ */
