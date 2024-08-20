
#ifndef FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_
#define FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_

#include <map>
#include <string>
#include <vector>
#include "math/include/histogram.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/macrostate.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class MonteCarlo;
class FlatHistogram;
class Clones;

/**
  Perform reweighting and ensemble averages using macrostate distributions
 */
class Ensemble {
 public:
  /// Store the original conjugate, macrostate and distribution.
  Ensemble(const Histogram& macrostates,
           const LnProbability& ln_prob,
           const double conjugate = 0.) { init_(macrostates, ln_prob); }

  /// Same as above, but taken from FlatHistogram.
  explicit Ensemble(const FlatHistogram& flat_hist);

  /// Same as above, but taken from spliced Clones.
  explicit Ensemble(const Clones& clones);

  /// Return the stored macrostate distribution from the original simulation.
  const LnProbability& ln_prob_original() const { return ln_prob_original_; }

  /// Return the LnProbability (that may have been reweighted).
  const LnProbability& ln_prob() const { return ln_prob_; }

  // HWH used for extrapoltaion.. deprecate?
  // Set the ln_prob.
  void set_ln_prob(const LnProbability& ln_prob) { ln_prob_ = ln_prob; }

  /// Return the stored Histogram representing the macrostates of the original
  /// simulation.
  const Histogram& macrostates() const { return macrostates_; }

  /// Determine min and max macrostate indices for a given phase
  /// Return 1 if no phase boundary, or 0 for success.
  void phase_boundary(const int phase, int * min, int * max) const;

  /// Return true if phase boundary is found.
  bool is_phase_boundary() const;

  /// Return the conjugate thermodynamic variable of the macrostate from the
  /// original simulation.
  double original_conjugate() const { return original_conjugate_; }

  /**
    Store and return a reweighted macrostate distribution due to a change
    in the conjugate variable.
    For example, if reweighting in number
    of particles using the grand canonical ensemble, the change in the
    thermodynamic conjugate is \f$\Delta(\beta\mu)\f$.
    If reweighting in potential energy in the microcanonical, the change
    in the thermodynamic conjugate is \f$\Delta(-\beta)\f$.
   */
  const LnProbability& reweight(const double delta_conjugate);

  /// Return the conjugate variable corresponding to the current (reweighted)
  /// macrostate distribution.
  double conjugate() const { return original_conjugate() + delta_conjugate_; }

  /// Return the ensemble averaged property from a list of properties averaged
  /// at each macrostate.
  double average(const std::vector<double>& macrostate_averages,
                 /// Select phase by order of macrostate.
                 /// Assumes default method of dividing phase boundary.
                 const int phase = 0) const;

  /// Return the average macrostate
  double average_macrostate(
    /// Select phase by order of macrostate.
    /// Assumes default method of dividing phase boundary.
    const int phase = 0) const;

 protected:
  void init_(const Histogram& macrostates,
             const LnProbability& ln_prob);
  double original_conjugate_ = 0.;
  LnProbability ln_prob_original_;

 private:
  Histogram macrostates_;
  LnProbability ln_prob_;
  double delta_conjugate_ = 0.;  // record reweight change in conjugate
};

/**
  Grand canonical ensemble currently implemented for single component systems.
 */
class GrandCanonicalEnsemble : public Ensemble {
 public:
  /// Store \f$\beta\f$, \f$\mu\f$ and the original ln_prob_.
  GrandCanonicalEnsemble(const Histogram& macrostates,
    const LnProbability& ln_prob,
    const double conjugate = 0.);

  // Same as above, but taken from FlatHistogram.
  GrandCanonicalEnsemble(const FlatHistogram& flat_histogram,
    const double beta_mu);

  // Same as above, but taken from spliced Clones.
  explicit GrandCanonicalEnsemble(const Clones& clones);

  /// Return the conjugate to the macrostate, \f$\beta\mu\f$.
  double beta_mu() const { return conjugate(); }

  /// Return \f$\beta PV\f$, the pressure times volume and inverse temperature.
  /// See https://doi.org/10.1063/1.2064628
  double betaPV(
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const;
};

class ExtrapolateBetaGCE : public GrandCanonicalEnsemble {
 public:
  /**
    Extrapolate the original simulation data to a different \f$\beta\f$.
    See: https://doi.org/10.1063/1.4975331
    Assumes a single component.
    In this implementation, extrapolation should be done before reweighting
    because this function changes ln_prob_original, not ln_prob.

    args:
    - beta_new: \f$\beta\f$ to extrapolate toward.
    - beta_original: original \f$\beta\f$ of the simulation.
    - order: truncate Taylor series at this order (default: 2).
   */
  ExtrapolateBetaGCE(const MonteCarlo& monte_carlo,
                     const FlatHistogram& flat_histogram,
                     argtype args = argtype());

  /// Same as above, but stitched together from parallel Clones.
  ExtrapolateBetaGCE(const Clones& clones,
    argtype args = argtype());

  /// Return the new energies based on extrapolation.
  const std::vector<double> energy() const { return energy_; }

 private:
  std::vector<double> energy_;

  void extrapolateBetaGCE_(
    /// Canonical ensemble average moments of extensive energy.
    /// The first index is the order of the moment.
    /// The second index is the macrostate.
    /// The first moment (corresponding to U) is modified to be the new energy
    /// at the new, extrapolated value of \f$\beta\f$.
    const std::vector<std::vector<double> >& energy_moments,
    argtype args = argtype());
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_
