#include "MultiLevelMonteCarlo.hpp"
#include "pnl/pnl_mathtools.h"
#include <iostream> 

MultiLevelMonteCarlo::MultiLevelMonteCarlo(Option *opt, Model *mod, PnlRng *rng)
{
  m_opt = opt;
  m_mod = mod;
  m_rng = rng;
}

long long MultiLevelMonteCarlo::nSamples(int level, int m, int L)
{
  return pow(m, 2*L-level)*L;
}


void MultiLevelMonteCarlo::collapse(PnlMat *Gcrude, PnlMat *Gfine, int m)
{
  pnl_mat_resize(Gcrude, Gfine->m/m, Gfine->n);
  for (int j = 0; j<Gcrude->n; j++){
    for (int i = 0; i<Gcrude->m; i++){
      pnl_mat_set(Gcrude, i, j, MGET(Gfine, m*i,j));
    }
  }
}

void MultiLevelMonteCarlo::run(double &prix, double &std_dev, int m, int L)
{
  PnlVect *pathCrude = pnl_vect_new();
  PnlVect *pathFine = pnl_vect_new();
  PnlMat *Gcrude = pnl_mat_new();
  PnlMat *Gfine = pnl_mat_new();
  double prix0;
  double std_dev0;
  double price;
  PnlMat *G = pnl_mat_create(m , 2);
  PnlVect *path = pnl_vect_create(m);

  // Treat Level 0
  for (int i = 0; i<nSamples(0, m, L); i++){
    pnl_mat_rng_normal(G, m, 2, m_rng);
    m_mod->simul(path, m_opt->m_maturity, m, G);
    price = m_opt->payoff(path);
    prix += price ;
    std_dev += price*price;
  }
  prix0 = prix/nSamples(0, m, L);
  std_dev0 = std::sqrt((std_dev /nSamples(0, m, L) - prix*prix)/nSamples(0, m, L));
  // Loop on all the levels > 0
  for (int l = 1; l < L; l++) {
    long long Nl = nSamples(l, m, L);
    int ml = pnl_pow_i(m, l);
    double sum = 0.;

    for (int i = 0; i < Nl; i++) {
      pnl_mat_rng_normal(Gfine, ml, 2, m_rng);

      // Simulation of the fine model
      m_mod->simul(pathFine, m_opt->m_maturity, ml, Gfine);
      sum += m_opt->payoff(pathFine);

      // Simulation of the crude model
      collapse(Gcrude, Gfine, m);

      m_mod->simul(pathCrude, m_opt->m_maturity, ml / m, Gcrude);

      sum -= m_opt->payoff(pathCrude);
    }
    sum /= Nl;
    prix += sum;
    std_dev += sum *sum;
  }
  prix = prix/L;
  std_dev = prix /L;
  pnl_mat_free(&Gcrude);
  pnl_mat_free(&Gfine);
  pnl_vect_free(&pathCrude);
  pnl_vect_free(&pathFine);
}
