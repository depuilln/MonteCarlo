
#include "Option.hpp"

PnlVect* brownien(PnlVect* normal, Option* option){
    double ecart = sqrt(option->T/option->J);
    PnlVect* gPrime = pnl_vect_create (option->J+1);
    pnl_vect_set (gPrime, 0, 0.0);
    for (int i = 1; i < option->J+1; i++){
        pnl_vect_set (gPrime, i, pnl_vect_get(gPrime, i-1)+ecart*pnl_vect_get(normal,i-1));
    }
    return gPrime;
}

void inversionNormal(PnlVect* normal){
    for(int i =0; i<normal->size;i++){
        pnl_vect_set(normal, i, -pnl_vect_get(normal, i));
    }
}

void monteCarlo_antithetique(double& price, double& std_dev, double M, Option* option){
    double somme1 = 0;
    double somme2 = 0;
    double cov = 0;
    double var1 = 0;
    double var2 = 0;

    for (int i = 0; i < M; i++) {
        PnlVect* g = pnl_vect_new();
        pnl_vect_rng_normal(g, option->J, option->rng);
        PnlVect* sousJacent1 = option->sousJacent(brownien(g, option));
        inversionNormal(g);
        PnlVect* sousJacent2 = option->sousJacent(brownien(g, option));


        somme1 += option->payOff(sousJacent1)* exp(-option->r*option->T);
        somme2 += option->payOff(sousJacent2)* exp(-option->r*option->T);

        var1 += SQR(option->payOff(sousJacent1)* exp(-option->r*option->T));
        var2 += SQR(option->payOff(sousJacent2)* exp(-option->r*option->T));

        cov += option->payOff(sousJacent1)* exp(-option->r*option->T)*option->payOff(sousJacent2)* exp(-option->r*option->T);
        
        free(g);
        free(sousJacent2);
        free(sousJacent1);
    }
    double somme = (somme1 + somme2) / 2/M;
    cov = cov - somme1/M * somme2/M;
    var1 = var1 / M - SQR(somme1/M);
    var2 = var2 / M - SQR(somme2/M);
    price = somme;
    std_dev = sqrt((var1 + var2 + 2*cov)/2/M);

}

void monteCarlo(double& price, double& std_dev, double M, Option* option){
    PnlVect* gPrime;
    double somme = 0;
    double estim = 0;
    for (int i = 0; i < M; i++) {
        PnlVect* g = pnl_vect_new();
        pnl_vect_rng_normal(g, option->J, option->rng);
        gPrime = brownien(g, option);

        PnlVect* sousJacent_ = option->sousJacent(gPrime);
        somme += option->payOff(sousJacent_)* exp(-option->r*option->T);
        estim += SQR(option->payOff(sousJacent_)* exp(-option->r*option->T));

        free(g);
        free(sousJacent_);
    }
    somme = somme/M;
    estim = estim/M - SQR(somme);

    price = somme;
    std_dev = sqrt(estim);
    free(gPrime);
}

void monteCarlo_controle(double& price, double& std_dev, double M, Option* option){
    double somme = 0;
    double cov = 0;
    double var1 = 0;
    double var2 = 0;

    for (int i = 0; i < M; i++) {
        PnlVect* g = pnl_vect_new();
        pnl_vect_rng_normal(g, option->J, option->rng);
        PnlVect* sousJacent_ = option->sousJacent(brownien(g, option));

        controle += pnl_vect_get(sousJacent_, option->J)-option->s0*exp(option->r*option->T);
        somme += option->payOff(sousJacent_)* exp(-option->r*option->T)-pnl_vect_get(sousJacent_, option->J)-option->s0*exp(option->r*option->T);

        somme +=
        var1 += SQR(option->payOff(sousJacent1)* exp(-option->r*option->T));
        var2 += SQR(option->payOff(sousJacent2)* exp(-option->r*option->T));

        cov += option->payOff(sousJacent1)* exp(-option->r*option->T)*option->payOff(sousJacent2)* exp(-option->r*option->T);
        
        free(g);
        free(sousJacent2);
        free(sousJacent1);
    }
    double varContr = somme/M - c * ()
    cov = cov - somme1/M * somme2/M;
    var1 = var1 / M - SQR(somme1/M);
    var2 = var2 / M - SQR(somme2/M);
    price = somme;
    std_dev = sqrt((var1 + var2 + 2*cov)/2/M);

}

int
main(int argc, char const* argv[]){
    Option* option = new Option(100.0, 2.0, 0.25, 0.02, 0.95, 24.0);
    double price;
    double std;
    monteCarlo(price, std, 5000.0, option);
    std::cout << "Monte Carlo : " << price << " : " << std << "\n";
    monteCarlo_antithetique(price, std, 5000.0, option);
    std::cout << "Antithetique : " << price << " : " << std << "\n";
    free(option);
}
