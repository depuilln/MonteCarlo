// vim: set sw=4 ts=4 sts=4:

#include "pnl/pnl_random.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_specfun.h"
#include "pnl/pnl_finance.h"


#include <iostream>
#include <cmath>
#include <math.h>

void brownien(PnlVect* g, PnlVect* gPrime, double N, double T, PnlRng* rng){
    double ecart = sqrt(T/N);
    pnl_vect_set (gPrime, 0, 0.0);
    for (int i = 1; i < N+1; i++){
        pnl_vect_set (gPrime, i, pnl_vect_get(gPrime, i-1)+ecart*pnl_vect_get(g,i-1));
    }
}

double payoff(double S0, double sigma, double r, double T, double K, double g){
    return std::max(S0 * std::exp((r-0.5*SQR(sigma))*T + sigma* std::sqrt(T)*g)-K, 0.);
}

double price(double S0, double sigma, double r, double T, double K, double g){
    return std::exp(-r*T) * payoff(S0, sigma, r, T, K, g);
}

void df(double &deltaDF, double S0, double sigma, double r, double T, double K, int N, PnlRng *rng)
{

    double eps = 1 / pow((double)N, 0.2);
    double mean1 = 0.;
    double mean2 = 0.;
    for (size_t i = 0; i < N; i++) {
        double g = pnl_rng_normal(rng);
        mean1 += price(S0+eps, sigma, r, T, K, g);
        mean2 += price(S0-eps, sigma, r, T, K, g);
    }
    deltaDF = (mean1 - mean2) / 2 /eps/N;
}

void vr(double &deltaVR, double S0, double sigma, double r, double T, double K, int N, PnlRng *rng)
{
    double mean1 = 0.;
    for (size_t i = 0; i < N; i++) {
        double g = pnl_rng_normal(rng);
        mean1 += price(S0, sigma, r, T, K, g)*g/sigma/S0/std::sqrt(T);
    }
    deltaVR = (mean1)/N;
}

int main()
{
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    double spot = 100.;
    double sigma = 0.2;
    double r = 0.05;
    double T = 2.;
    double K = 110.;
    double g = pnl_rng_normal (rng);
    double N = 50000;


    double price, delta;
    pnl_cf_call_bs (spot, K, T, r, 0., sigma, &price, &delta);

    double deltaDF;
    df(deltaDF, spot, sigma, r, T, K, N, rng);

    double deltaVR;
    vr(deltaVR, spot, sigma, r, T, K, N, rng);

    std::cout << delta << " : Delta\n";
    std::cout << deltaDF << " : DeltaDF\n";
    std::cout << deltaVR << " : DeltaVR\n";
}
