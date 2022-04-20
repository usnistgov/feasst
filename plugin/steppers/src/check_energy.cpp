#include <cmath>
#include "utils/include/serialize.h"
#include "steppers/include/check_energy.h"

namespace feasst {

class MapCheckEnergy {
 public:
  MapCheckEnergy() {
    CheckEnergy().deserialize_map()["CheckEnergy"] = MakeCheckEnergy();
  }
};

static MapCheckEnergy mapper_energy_check_ = MapCheckEnergy();

CheckEnergy::CheckEnergy(argtype * args) : ModifyUpdateOnly(args) {
  tolerance_ = dble("tolerance", args, 1e-10);
  check_ = MakeCheck();
}
CheckEnergy::CheckEnergy(argtype args) : CheckEnergy(&args) {
  check_all_used(args);
}

void CheckEnergy::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  check_->update(*criteria, *system, *trial_factory);
  DEBUG("computing unoptimized energy for check");

  const double energy = system->unoptimized_energy();
  const double current_energy = criteria->current_energy();
  DEBUG("energy:" << energy << " "
     << "current_energy: " << current_energy << " "
     << "diff: " << energy - current_energy
  );
  accumulator_.accumulate(energy - current_energy);

  // loop over each profile and perform the energy check
  const std::vector<double>& energy_profile = system->unoptimized().stored_energy_profile();
  const std::vector<double>& current_energy_profile = criteria->current_energy_profile();
  DEBUG("energy_profile " << feasst_str(energy_profile));
  DEBUG("current_energy_profile " << feasst_str(current_energy_profile));
  for (int i = 0; i < static_cast<int>(energy_profile.size()); ++i) {
    ASSERT(std::abs(energy_profile[i] - current_energy_profile[i]) < tolerance_,
      MAX_PRECISION <<
      "Energy check failure. There is a problem with the potentials. " <<
      "The unoptimized energy of the potential " << i << " was computed as " <<
      energy_profile[i] << " but the running energy from criteria " <<
      "(the accumulation of a change in energy over a series of steps) is " <<
      current_energy_profile[i] <<
      ". The difference(" << std::abs(energy_profile[i] - current_energy_profile[i]) << ") is " <<
      "greater than the tolerance(" << tolerance_ << "). ");
  }
  criteria->set_current_energy_profile(energy_profile);

  // perform same energy check for entire system
  ASSERT(std::abs(energy - current_energy) < tolerance_,
    MAX_PRECISION <<
    "Energy check failure. There is a problem with the potentials. " <<
    "The unoptimized energy of the entire configuration was computed as " <<
    energy << " but the running energy from criteria " <<
    "(the accumulation of a change in energy over a series of steps) is " <<
    current_energy <<
    ". The difference(" << std::abs(energy - current_energy) << ") is " <<
    "greater than the tolerance(" << tolerance_ << "). "
    << system->unoptimized().str());
  criteria->set_current_energy(energy);

  // loop over all queryable maps and check those as well.
}

void CheckEnergy::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(715, ostr);
  feasst_serialize(tolerance_, ostr);
  feasst_serialize_fstdr(check_, ostr);
}

CheckEnergy::CheckEnergy(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 715, "version mismatch: " << version);
  feasst_deserialize(&tolerance_, istr);
  // feasst_deserialize_fstdr(modify->check_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      check_ = check_->deserialize(istr);
    }
  }
}

}  // namespace feasst
