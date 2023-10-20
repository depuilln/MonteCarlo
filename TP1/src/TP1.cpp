#include "pnl/pnl_random.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_specfun.h"
#include <iostream>
#include <cmath>
#include <math.h>

void brownien(PnlVect* g, PnlVect* gPrime, double N, double T, PnlRng* rng){
    pnl_vect_rng_normal(g, N, rng);
    double ecart = sqrt(T/N);
    pnl_vect_set (gPrime, 0, 0.0);
    for (int i = 1; i < N+1; i++){
        pnl_vect_set (gPrime, i, pnl_vect_get(gPrime, i-1)+ecart*pnl_vect_get(g,i-1));
    }
}

double estim_plus_precis(double a, double N, double T){
    PnlRng* rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    PnlVect* g = pnl_vect_new();
    PnlVect* gPrime = pnl_vect_create (N+1);
    double somme = 0;
    bool condition;

    for (int i = 0; i < 5000; i++) {
        brownien(g, gPrime, N, T, rng);
        int h = 0;
        for (int j = 0; j<N; j++){

            condition = a > 0 ? pnl_vect_get(gPrime,h)<a : pnl_vect_get(gPrime,h)>a;
            if (!condition){
                double t = h;
                double s = h-T;
                for (int k = 1; k<=10; k++){
                    double u = (h-T) + k*T/10.0;

                    double nouveau = (t-u)/(t-s)*pnl_vect_get(gPrime,h-T)+(u-s)/(t-s)*pnl_vect_get(gPrime,h)+ sqrt((t-u)*(u-s)/(t-s)) * pnl_rng_normal (rng);
                    if (nouveau>a){
                        h = u;
                    }
                }
                break;
            }
            else {
                h+=1;
                condition = a > 0 ? pnl_vect_get(gPrime,h)<a : pnl_vect_get(gPrime,h)>a;
            }
        }
        somme += h*T/N;
    }
    return somme/ 5000.0;

    free(g);
    free(gPrime);
    free(rng);
}

double estim(double a, double N, double T){
    PnlRng* rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    PnlVect* g = pnl_vect_new();
    PnlVect* gPrime = pnl_vect_create (N+1);
    double somme = 0;
    for (int i = 0; i < 2000; i++) {
        brownien(g, gPrime, N, T, rng);

        int h = 0;
        bool condition = a > 0 ? pnl_vect_get(gPrime,h)<a : pnl_vect_get(gPrime,h)>a;
        while(condition and h<N){
            h+=1;
            condition = a > 0 ? pnl_vect_get(gPrime,h)<a : pnl_vect_get(gPrime,h)>a;
        }
        somme += h*T/N;
    }
    return somme/ 2000.0;

    free(g);
    free(gPrime);
    free(rng);
}

double valeur_theorique(double a, double T){
    double gamma1 = pnl_sf_gamma_inc(0.5, SQR(a)/2.0/T);
    double gamma2 = pnl_sf_gamma_inc(-0.5, SQR(a)/2.0/T);
    return T*(1-gamma1/sqrt(M_PI))+gamma2*SQR(a)/2.0/sqrt(M_PI);
}

int
main(int argc, char const* argv[]){
    
    double a = 2.0;
    // for (int i = 1; i<6; i++){
    //     std::cout << "T = " << i << " : " << estim_plus_precis(a, 500*i, i) << ", " << valeur_theorique(a, i) << "\n";
    // }

    for (int i = 1; i<50; i++){
        std::cout << (estim_plus_precis(a, 500*i, 5)-valeur_theorique(a, 5))*sqrt(500*i) << "\n";
    }
}
