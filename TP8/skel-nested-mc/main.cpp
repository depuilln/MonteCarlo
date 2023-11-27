#include <ctime>
#include <iostream>

#include "AsianOption.hpp"
#include "BlackScholes.hpp"
#include "MonteCarlo.hpp"

int main(int argc, char **argv)
{
  BlackScholesModel bs = BlackScholesModel(0.03, 100., 0.3);
  AsianOption asianOption = AsianOption(110, 2);
  PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
  pnl_rng_sseed(rng, std::time(NULL));
  MonteCarlo mc = MonteCarlo(&asianOption, &bs, rng);

  double prix, prix_imbr, std_dev, std_dev_imbr;
  int nTimeSteps = 50;
  long long nSamples = 1E5;
  long mSamples = 20;
  mc.run(prix, std_dev, nSamples, nTimeSteps);
  std::cout << "Price: " << prix << "\n";
  std::cout << "CI width: " << std_dev * 1.96 * 2 << "\n";

  mc.runNested(prix_imbr, std_dev_imbr, nSamples, nTimeSteps, 0.8, 1.);
  std::cout << "Price Imbr: " << prix_imbr << "\n";
  std::cout << "CI width Imbr: " << std_dev_imbr * 1.96 * 2 << "\n";

  pnl_rng_free(&rng);
  return 0;
}