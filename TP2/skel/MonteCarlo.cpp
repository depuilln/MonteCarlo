// vim: set sw=4 ts=4 sts=4:
#include "pnl/pnl_random.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_specfun.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "MonteCarlo.hpp"

MonteCarlo::MonteCarlo(BSCall &product, size_t samples)
    : m_product(product), m_samples(samples)
{}



void MonteCarlo::mc(double &prix, double &stddev, PnlRng *rng)
{
    double sumX = 0.;
    double var = 0.;
    for (size_t i = 0; i < m_samples; i++) {
        double g = pnl_rng_normal(rng);
        double flow = m_product.payoff(g);
        sum += flow;
        var += flow * flow;
    }
    prix = sum / m_samples;
    var = var / m_samples - prix * prix;
    stddev = std::sqrt(var / m_samples);
}

void MonteCarlo::is(PnlVect *lambda, double gamma, int n, PnlRng *rng){
    double beta = 0.75;
    PnlVect* g = pnl_vect_new();
    pnl_vect_rng_normal(g, n, rng);
    double alpha = 0;
    double theta = 0;
    for (int i = 1; i<=n; i++){
        theta = theta - (gamma/pow(i+1, beta))*this->m_product.d_weight_plus(theta, pnl_vect_get(g,i-1)); 
        if (SQR(theta)> log(alpha+1)){
            theta = 0;
            alpha++;
        }
        pnl_vect_set(lambda, i-1, theta);
    }
}

void MonteCarlo::mcis(double &prix, double &stddev, double lambda, PnlRng *rng){
    double sum = 0.;
    double var = 0.;
    for (size_t i = 0; i < m_samples; i++) {
        double g = pnl_rng_normal(rng);
        double flow = m_product.payoff(g, lambda)*m_product.weight_minus(g, lambda);
        sum += flow;
        var += flow * flow;
    }
    prix = sum / m_samples;
    var = var / m_samples - prix * prix;
    stddev = std::sqrt(var / m_samples);
}

