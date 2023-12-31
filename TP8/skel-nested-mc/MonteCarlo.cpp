#include <cmath>
#include "MonteCarlo.hpp"
#include "pnl/pnl_mathtools.h"

MonteCarlo::MonteCarlo(Option *opt, BlackScholesModel *mod, PnlRng *rng)
{
  m_opt = opt;
  m_mod = mod;
  m_rng = rng;
}

void MonteCarlo::run(double &prix, double &std_dev, long long nSamples, int nTimeSteps)
{
  PnlVect *G = pnl_vect_new();
  PnlVect *path = pnl_vect_new();
  double sum = 0.;
  double sum_sq = 0.;
  for (int n = 0; n < nSamples; n++) {
    pnl_vect_rng_normal(G, nTimeSteps, m_rng);
    m_mod->simul(path, m_opt->m_maturity, nTimeSteps, G);
    double payoff = m_opt->payoff(path);
    sum += payoff;
    sum_sq += payoff * payoff;
  }
  prix = std::exp(- m_mod->m_r * m_opt->m_maturity) * sum / nSamples;
  std_dev = std::sqrt((std::exp(- 2 * m_mod->m_r * m_opt->m_maturity) * sum_sq / nSamples - prix * prix) / nSamples);
  pnl_vect_free(&G);
  pnl_vect_free(&path);
}

void MonteCarlo::runNested(double &prix, double &std_dev, long long nSamples, int nTimeSteps, double alpha, double t)
{
  PnlVect *G = pnl_vect_new();
  PnlVect *path = pnl_vect_new();
  double sum = 0.;
  double sum_sq = 0.;
  double portT, portT_sq, port0;
  for (int n = 0; n < nSamples; n++) {
    run(port0, portT_sq, nSamples, nTimeSteps);
    runT(portT, portT_sq, nSamples, nTimeSteps);
    double payoff = MAX(portT - alpha * port0, 0);
    sum += payoff;
    sum_sq += payoff * payoff;
  }
  prix = std::exp(- m_mod->m_r * m_opt->m_maturity) * sum / nSamples;
  std_dev = std::sqrt((std::exp(- 2 * m_mod->m_r * m_opt->m_maturity) * sum_sq / nSamples - prix * prix) / nSamples);
  pnl_vect_free(&G);
  pnl_vect_free(&path);
}

void MonteCarlo::runT(double &prixPorto, double &std_devPorto, long long mSamples, int nTimeSteps, double t)
{
  PnlVect *G = pnl_vect_new();
  PnlVect *path = pnl_vect_new();
  double sum = 0.;
  double sum_sq = 0.;
  for (int n = 0; n < mSamples; n++){
    pnl_vect_rng_normal(G, nTimeSteps, m_rng);
    m_mod->simul(path, m_opt->m_maturity, nTimeSteps, G);
    double payoff = m_opt->payoff(path);
    sum += payoff;
    sum_sq += payoff * payoff;
  }
  prixPorto = std::exp(- m_mod->m_r * (m_opt->m_maturity-t))*sum/mSamples;
  std_devPorto = std::sqrt((std::exp(- 2 * m_mod->m_r * (m_opt->m_maturity-t)) * sum_sq / nSamples - prix * prix) / nSamples);
  pnl_vect_free(&G);
  pnl_vect_free(&path);
}