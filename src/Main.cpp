#include "../include/Utils.h"
#include "../include/BackgroundCosmology.h"
#include "../include/RecombinationHistory.h"
#include "../include/Perturbations.h"
#include "../include/PowerSpectrum.h"
#include "../include/SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB0     = 0.05;
  double OmegaCDM0   = 0.267;
  double OmegaK0     = 0.0;
  double Neff        = 3.046;
  double TCMB0       = 2.7255;

  // Recombination parameters
  double Yp          = 0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB0);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("data/cosmology.txt");

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  // Utils::StartTiming("Supernova Fitting");
  mcmc_fit_to_supernova_data("data/supernovadata.txt", "data/results_supernovafitting.txt");
  // // End timing for supernova fitting
  // Utils::EndTiming("Supernova Fitting");

  // Remove when module is completed
  Utils::EndTiming("Everything");
  return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("data/recombination.txt");
  
  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "data/perturbations_k0.01.txt");
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("data/cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
