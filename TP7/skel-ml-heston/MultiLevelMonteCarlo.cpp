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
  double price = 0.;
  PnlMat *G = pnl_mat_create(m , 2);
  PnlVect *path = pnl_vect_create(m);

  // Treat Level 0
  for (int i = 0; i<nSamples(0, 1, L); i++){
    pnl_mat_rng_normal(G, m, 2, m_rng);
    m_mod->simul(path, m_opt->m_maturity, m, G);
    price = m_opt->payoff(path);
    prix += price ;
    std_dev += price*price;
  }
  prix /= nSamples(0, 1, L);
  std_dev /= nSamples(0, 1, L);
  std_dev = std::sqrt((std_dev - prix*prix)/nSamples(0, 1, L));
  // Loop on all the levels > 0
  for (int l = 1; l < L; l++) {
    long long Nl = nSamples(l, m, L);
    int ml = pnl_pow_i(m, l);
    double sum = 0.;
    double var = 0.;
    for (int i = 0; i < Nl; i++) {
      pnl_mat_rng_normal(Gfine, ml, 2, m_rng);

      // Simulation of the fine model
      m_mod->simul(pathFine, m_opt->m_maturity, ml, Gfine);
      sum += m_opt->payoff(pathFine);

      // Simulation of the crude model
      collapse(Gcrude, Gfine, m);

      m_mod->simul(pathCrude, m_opt->m_maturity, ml / m, Gcrude);

      sum -= m_opt->payoff(pathFine);
      var += SQR(m_opt->payoff(pathFine)-m_opt->payoff(pathFine));
    }
    sum /= Nl;
    var /= Nl;
    prix += sum;
    std_dev += std::sqrt((var - sum*sum)/Nl);
  }
  pnl_mat_free(&Gcrude);
  pnl_mat_free(&Gfine);
  pnl_vect_free(&pathCrude);
  pnl_vect_free(&pathFine);
  pnl_mat_free(&G);
  pnl_vect_free(&path);
}
