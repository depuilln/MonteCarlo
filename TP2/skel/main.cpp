// vim: set sw=4 ts=4 sts=4:

#include <iostream>
#include <ctime>
#include "MonteCarlo.hpp"
#include "BSCall.hpp"
#include <fstream>
int main()
{
    std::ofstream output("frere.ods");

    double prix, stddev;
    BSCall product(2., 0.2, 0.03, 100., 120.);
    MonteCarlo pricer(product, 50000);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, std::time(NULL));
    pricer.mc(prix, stddev, rng);
    std::cout << "prix : " << prix << " (IC = " << stddev * 1.96 << ")\n";

    int n = 200;
    double gamma = 50;
    PnlVect *theta = pnl_vect_create(n);
    pricer.is(theta, gamma, n, rng);
    for (int i = 0; i<n; i++){
        output << pnl_vect_get(theta, i)<<"\n";
    }
    pricer.mcis(prix, stddev, pnl_vect_get(theta, n-1), rng);
    std::cout << "prix : " << prix << " (IC = " << stddev * 1.96 << ")\n";
    output.close();
    pnl_rng_free(&rng);
    exit(0);
}
