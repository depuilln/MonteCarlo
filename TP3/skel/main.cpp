// vim: set sw=4 ts=4 sts=4:

#include "pnl/pnl_random.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_specfun.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include "MonteCarlo.hpp"
#include "BSBarrier.hpp"

void brownien(PnlVect* g, PnlVect* gPrime, double N, double T, PnlRng* rng){
    double ecart = sqrt(T/N);
    pnl_vect_set (gPrime, 0, 0.0);
    for (int i = 1; i < N+1; i++){
        pnl_vect_set (gPrime, i, pnl_vect_get(gPrime, i-1)+ecart*pnl_vect_get(g,i-1));
    }
}

double u(double lambda, PnlMat* W, PnlVect* fcarre, double m_dates){
    
    double somme1 = 0;
    double somme2 = 0;
    double somme3 = 0;
    double wi;
    for (int i = 0; i<fcarre->size; i++){
        wi = MGET(W, i, W->n-1);
        somme1 += wi*exp(-lambda * wi)*pnl_vect_get(fcarre, i);
        somme2 += exp(-lambda * wi)*pnl_vect_get(fcarre, i);
        somme3 += SQR(wi)*exp(-lambda * wi)*pnl_vect_get(fcarre, i);
    }
    double uPrime = lambda - somme1/somme2;
    double uSec = 1 + (somme3*somme2-SQR(somme1))/SQR(somme2);
    return uPrime/uSec;
}

double lambdaCV(int p, int n, PnlMat* W, PnlVect* fcarre, double m_dates){
    double lambda = 0;
    for (int i = 0; i<p; i++){
        lambda -= u(lambda, W, fcarre, m_dates);
    }
    return lambda;
}

int main()
{
    double m_maturity = 2.;
    double m_volatility = 0.2;
    double m_interest_rate = 0.05;
    double m_spot = 100.;
    double m_strike = 110.;
    double m_barrier = 80.;
    int m_dates = 24;
    double prix, stddev;
    BSBarrier product(m_maturity, m_volatility, m_interest_rate, m_spot, m_strike, m_barrier, m_dates);
    MonteCarlo pricer(product, 50000);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    pricer.mc(prix, stddev, rng);
    std::cout << "prix : " << prix << " (half-IC = " << stddev * 1.96 << ")\n";

    int N = 100;
    int n = 500;

    PnlMat* W = pnl_mat_create (n, m_dates+1);
    PnlVect* brow = pnl_vect_create (m_dates+1);
    PnlVect* g; 
    PnlVect* fcarre = pnl_vect_create(n);
    for (int i = 0; i<n; i++){
        g = pnl_vect_create(m_dates);
        pnl_vect_rng_normal(g, m_dates, rng);
        brownien(g, brow, m_dates, m_maturity, rng);
        PnlVect* sousJacent = pnl_vect_create(m_dates+1);
        product.asset(sousJacent, brow);
        pnl_mat_set_row (W, brow, i);

        pnl_vect_set(fcarre, i, SQR(product.payoff(sousJacent)));
    }

    double lambda = lambdaCV(10, n, W, fcarre, m_dates);
    pricer.mcVar(prix, stddev, lambda, rng);
    std::cout << "prix : " << prix << " (half-IC = " << stddev * 1.96 << ")\n";

    free(fcarre);
    free(brow);
    free(g);
    free(W);
    pnl_rng_free(&rng);
}
