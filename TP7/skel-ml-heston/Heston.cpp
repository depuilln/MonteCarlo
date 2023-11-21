#include "Heston.hpp"
#include <iostream>
#include <cmath>
#include "pnl/pnl_vector.h"

HestonModel::HestonModel(double r, double spot, double initVol, double sigma, double kappa, double theta, double rho)
  : Model(2, r)
  , m_spot(spot)
  , m_initVol(initVol)
  , m_sigma(sigma)
  , m_kappa(kappa)
  , m_theta(theta)
  , m_rho(rho)
{}

void HestonModel::simul(PnlVect *path, double maturity, double nTimeSteps, PnlMat *G)
{
  pnl_vect_resize(path, nTimeSteps + 1);
  pnl_vect_set(path, 0, m_spot);
  double dt = maturity / nTimeSteps;
  double v = m_initVol;
  double value = 0;
  for (int i = 1; i< nTimeSteps +1; i++){
    v = v + m_kappa * (m_theta - v) * dt + m_sigma * std::sqrt(std::max(v,0.)) * MGET(G, i-1, 1)*std::sqrt(dt);
    value = GET(path, i-1);
    pnl_vect_set(path, i, value+m_r*value*dt+std::sqrt(v)*value*MGET(G, i-1,0)*std::sqrt(dt));
  }
  
}