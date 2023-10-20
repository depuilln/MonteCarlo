#pragma once
#include "pnl/pnl_random.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_specfun.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Classe Option abstraite
class Option
{
  public:
    PnlRng* rng;
    double T;        /// maturité
    double sigma;
    double r;
    double s0;
    double K;
    double J;

    Option(double s0, double T, double sigma, double r, double K, double J);
    /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d : N = Timestep number (nbTimeSteps_) et d : option size (size_)
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset.
     * @param[in] option est l'option duquel on calcule le pay off
     * @return payoff correspondant à la valeur du payoff
     */
    double payOff(PnlVect* sousJacent);

    PnlVect* sousJacent(PnlVect* brownien);

};
