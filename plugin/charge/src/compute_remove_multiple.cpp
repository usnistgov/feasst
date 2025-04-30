#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "charge/include/compute_remove_multiple.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

ComputeRemoveMultiple::ComputeRemoveMultiple(argtype args)
  : ComputeRemoveMultiple(&args) {
  feasst_check_all_used(args);
}
ComputeRemoveMultiple::ComputeRemoveMultiple(argtype * args)
  : TrialCompute(args) {
  class_name_ = "ComputeRemoveMultiple";
  shift_ = integer("shift", args, -1);
}

FEASST_MAPPER(ComputeRemoveMultiple,);

void ComputeRemoveMultiple::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeRemoveMultiple");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  if (!acceptance->reject()) {
    ASSERT(system->num_configurations() == 1, "not implemented for multiple configs");
    acceptance->set_energy_new(criteria->current_energy() - acceptance->energy_old());
    acceptance->set_energy_profile_new(criteria->current_energy_profile());
    acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old());
    acceptance->add_to_macrostate_shift(shift_);
    acceptance->set_macrostate_shift_type(-1);  // disable constraints with multi-particles
    DEBUG("deltaE " << MAX_PRECISION << -1.*acceptance->energy_old());
    const Configuration& config = system->configuration();
    const double volume = config.domain().volume();

    // initialize delta_
    const int num_ptypes = config.num_particle_types();
    if (static_cast<int>(delta_.size()) < num_ptypes) delta_.resize(num_ptypes);
    std::fill(delta_.begin(), delta_.end(), 0);

    for (const TrialStage * stage : *stages) {
      const TrialSelect& select = stage->trial_select();
      const int particle_index = select.mobile().particle_index(0);
      const int particle_type = config.select_particle(particle_index).type();
      DEBUG("lnmet " << acceptance->ln_metropolis_prob());
      DEBUG("volume " << volume << " selprob " << select.probability()
        << " betamu " << system->thermo_params().beta_mu(particle_type));
      const int num_pt = config.num_particles_of_type(particle_type);
      const double prob = 1./static_cast<double>(num_pt - delta_[particle_type]);
      ++delta_[particle_type];
      acceptance->add_to_ln_metropolis_prob(
        - std::log(volume*prob)
        //- std::log(volume*select.probability())
        - system->thermo_params().beta_mu(particle_type)
      );
    }
  }
}

std::shared_ptr<TrialCompute> ComputeRemoveMultiple::create(
    std::istream& istr) const {
  return std::make_shared<ComputeRemoveMultiple>(istr);
}

ComputeRemoveMultiple::ComputeRemoveMultiple(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeRemoveMultiple", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(294 == version, "mismatch version: " << version);
  feasst_deserialize(&shift_, istr);
}

void ComputeRemoveMultiple::serialize_compute_remove_multiple_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(294, ostr);
  feasst_serialize(shift_, ostr);
}

void ComputeRemoveMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_remove_multiple_(ostr);
}

}  // namespace feasst
