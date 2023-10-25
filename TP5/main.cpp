#include "pnl/pnl_random.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_specfun.h"
#include "pnl/pnl_finance.h"

#include "MonteCarlo.hpp"
#include "BSCall.hpp"
#include <iostream>
#include <cmath>
#include <math.h>



int main()
{
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    double prix, prixY, prixZ, stddev, stddevY, stddevZ;
    BSCall product(1., 0.2, 0.095, 100., 100., 24);
    MonteCarlo pricer(product, 50000);
    
    pricer.mc(prix, prixY, prixZ, stddev, stddevY, stddevZ, rng);
    std::cout << "prixX : " << prix << " ICX = " << stddev * 1.96 << "\n";
    std::cout << "prixY : " << prixY <<" ICY : " << stddevY * 1.96 << "\n";
    std::cout << "prixZ : " << prixZ <<" ICZ : " << stddevZ * 1.96 << "\n";



    // double price, delta;
    // pnl_cf_call_bs (spot, K, T, r, 0., sigma, &price, &delta);

    // double deltaDF;
    // df(deltaDF, spot, sigma, r, T, K, N, rng);

    // double deltaVR;
    // vr(deltaVR, spot, sigma, r, T, K, N, rng);

    // std::cout << delta << " : Delta\n";
    // std::cout << deltaDF << " : DeltaDF\n";
    // std::cout << deltaVR << " : DeltaVR\n";
}