
#ifndef FEASST_GIBBS_GIBBS_INITIALIZE_H_
#define FEASST_GIBBS_GIBBS_INITIALIZE_H_

#include <vector>
#include <memory>
#include "monte_carlo/include/modify.h"

namespace feasst {

class TrialAdd;
class TrialRemove;
class TrialVolume;

typedef std::map<std::string, std::string> argtype;

/**
  This class attempts to initialize a Gibbs simuilation.
  Here we refer to the low density as vapor and high density as liquid.
  The number of particles and volume in the liquid is N_L and V_L, while N_V and
  V_V are in the vapor.
  Because the densities of the vapor and liquid may be unknown at the start of
  a simulation, this class attempts to add/remove volume or particles from the
  lower density configuration to obtain two target conditions.

  The equations are as follows:

  \f$N = \rho_V V_V + \rho_L V_L\f$

  \f$V = V_V + V_L\f$

  This is two equations with four unknowns.
  Therefore, we should introduce two constraints.

  The first constraint is that a given fraction of the total number of particles
  should be in the vapor phase.
  This leads to an equation as follows:

  \f$\rho_V V_V = f N\f$,

  where f is the given fraction in this example.
  This sets the relative size of the two configurations.

  The second constraint is that the number of particles will remain fixed.
  In this case,

  \f$V_L = N(1-f)/\rho_L\f$

  \f$V_V = Nf/\rho_v\f$

  and the low density Configuration is adjusted yield the predicted total
  volume.

  When all constraints pass within their given tolerance, the Criteria is set
  to complete to allow for production simulations.
 */
class GibbsInitialize : public Modify {
 public:
  //@{
  /** @name Arguments
    - particle_type: type of particle to consider density (default: 0).
    - fraction_particles_low_density: fractional amount of particles in the low
      density phase (default: 0.15).
    - fraction_particles_low_density_tolerane: tolerance on the fractional
      amount of particles in the low density phase (default: 0.05).
      When all tolerances are achieved, the Criteria is set to complete.
    - updates_density_equil: number of updates before computing average density
      (default: 1e7).
    - updates_per_adjust: number of updates per adjustment (default: 2e7).
   */
  explicit GibbsInitialize(argtype args = argtype());
  explicit GibbsInitialize(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override;

  std::string write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("GibbsInitialize"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<GibbsInitialize>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<GibbsInitialize>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit GibbsInitialize(std::istream& istr);
  ~GibbsInitialize();

  //@}
 private:
  double frac_part_low_dens_;
  int low_dens_config_;
  int high_dens_config_;
  int particle_type_;
  int updates_since_adjust_ = 0;
  int updates_density_equil_;
  int updates_per_adjust_;
  double frac_part_low_dens_tol_;

  std::unique_ptr<Accumulator> low_dens_;
  std::unique_ptr<Accumulator> high_dens_;

  int dens_config_index_(const bool high, // high dens if true
    const System& system) const;
};

inline std::shared_ptr<GibbsInitialize> MakeGibbsInitialize(
    argtype args = argtype()) {
  return std::make_shared<GibbsInitialize>(args);
}

}  // namespace feasst

#endif  // FEASST_GIBBS_GIBBS_INITIALIZE_H_
