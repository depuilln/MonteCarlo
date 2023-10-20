#include "Option.hpp"

Option::Option(double s0, double T, double sigma, double r, double K, double J)
{
    this->T = T;
    this->sigma = sigma;
    this->r = r;
    this->s0 = s0;
    this->K = K;
    this->J = J;
    this->rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(this->rng, time(NULL));
}

double Option::payOff(PnlVect* sousJacent){
    return std::max(K*pnl_vect_max(sousJacent) - pnl_vect_get(sousJacent, J), 0.0);
}

PnlVect* Option::sousJacent(PnlVect* brownien){
    PnlVect* traj = pnl_vect_create (J+1);
    pnl_vect_set(traj, 0, s0);
    for (int i = 1; i<=J; i++){
        double value = pnl_vect_get(brownien, i);
        pnl_vect_set(traj, i, s0 * exp((r - SQR(sigma)/2)*i*T/J+sigma*value));
    }
    return traj;
}