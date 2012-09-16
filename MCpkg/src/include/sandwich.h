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
#include <scythestat/la.h>
#include <cmath>
#include <omp.h>

namespace sandwich {
    namespace {
       typedef unsigned int uint;
    }

    // forward declarations
    scythe::Matrix<> estimationFunction(uint& index, const scythe::Matrix<>& data, const scythe::Matrix<>& par);
    scythe::Matrix<> computeWeights(const double& stepsize, const uint& pos_w);

    inline scythe::Matrix<>
    bread(const scythe::Matrix<>& cov_m, const uint& nobs, const double& ssr) {

    	const uint nc_m = cov_m.cols();
    	scythe::Matrix<> bread_m(cov_m);
    	double sigma2 = ssr / (nobs  - nc_m);
    	bread_m = bread_m * sigma2 * nobs;

        return bread_m;
    }

    inline scythe::Matrix<>
    meat(const scythe::Matrix<>& cov_m, const scythe::Matrix<>& data, const scythe::Matrix<>& par) {

    	const uint nobs = data.rows();
    	const uint nreg = data.cols() - 1;
    	double bandwidth = std::pow((double)nobs, 0.2);
        double stepsize = 1.0 / (1.0 + bandwidth);
    	uint pos_w = 1 + std::floor(1.0 / stepsize);

    	scythe::Matrix<> weights = computeWeights(stepsize, pos_w);
    	scythe::Matrix<> estimationf(nobs, nreg);

        #pragma omp parallel for // TODO: Check if private(data, par) is needed or an ordered instruction
    	for(uint i = 0; i < nobs; ++i) {
    		estimationf(i, scythe::_) = estimationFunction(i, data, par);
    	}
    	scythe::Matrix<> meat_m(nreg, nreg);
    	meat_m = t(estimationf) * estimationf * 0.5;

    	if(weights.rows() > 1) {
    		for(uint i = 0; i < weights.cols(); ++i) {
    			meat_m += (t(estimationf(0, 0, nobs - i, nreg)) * estimationf(i, 0, nobs - i, nreg));
    		}
    	}

    	double dim = 1.0/nobs;
    	meat_m += t(meat_m);
    	meat_m *= dim;

    	return meat_m;
    }

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

    inline scythe::Matrix<>
    estimationFunction(uint& index, const scythe::Matrix<>& data, const scythe::Matrix<>& par) {

    	const uint nc_d = data.cols();
        scythe::Matrix<> tmp_dep = data(index, 0, index, 0);
        scythe::Matrix<> tmp_hat = data(index, 1, index, 1) * (par(0,0) + par(1,0)) + data(index, 2, index, 2) * (par(0,0) + par(2, 0) * par(1,0));
        tmp_hat -= tmp_dep;
        tmp_hat *= tmp_hat * (-1);

        return data(index, 1, index, nc_d) * tmp_hat;
    }
}


#endif /* SANDWICH_H_ */
