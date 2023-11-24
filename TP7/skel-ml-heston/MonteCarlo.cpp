#include "MonteCarlo.hpp"
#include <iostream>
#include <cmath>

MonteCarlo::MonteCarlo(Option *opt, Model *mod, PnlRng *rng)
{
  m_opt = opt;
  m_mod = mod;
  m_rng = rng;
}

void MonteCarlo::run(double &prix, double &std_dev, long long nSamples, int nTimeSteps)
{
  PnlMat *G = pnl_mat_create(nTimeSteps , 2);
  PnlVect *path = pnl_vect_create(nTimeSteps);
  double price;
  for (int i = 0; i<nSamples; i++){
    pnl_mat_rng_normal(G, nTimeSteps, 2, m_rng);
    m_mod->simul(path, m_opt->m_maturity, nTimeSteps, G);

    price = m_opt->payoff(path);
    prix += price;
    std_dev += price*price;
  }
  pnl_mat_free(&G);
  pnl_vect_free(&path);
  prix = prix/nSamples;
  std_dev = std::sqrt((std_dev /nSamples - prix*prix)/nSamples);
}

void MonteCarlo::mse(double &err, long long nSamples, int nTimeSteps, double prixExact){
  double prix = 0.;
  double std_dev = 0.;
  for (int i = 0; i<20; i++){
    run(prix, std_dev, nSamples, nTimeSteps);
    err += (prix-prixExact)*(prix-prixExact);
  }
  err /=20;
}