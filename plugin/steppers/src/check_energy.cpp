#include <cmath>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "steppers/include/check_energy.h"

namespace feasst {

FEASST_MAPPER(CheckEnergy,);

CheckEnergy::CheckEnergy(argtype * args) : ModifyUpdateOnly(args) {
  tolerance_ = dble("tolerance", args, 1e-10);
  decimal_places_ = integer("decimal_places", args, -1);
  check_ = MakeCheck();
}
CheckEnergy::CheckEnergy(argtype args) : CheckEnergy(&args) {
  feasst_check_all_used(args);
}

bool CheckEnergy::is_within_tolerance_(const double u1, const double u2) const {
  bool pass = true;
  DEBUG(MAX_PRECISION << "u1 " << u1 << " u2 " << u2);
  if (decimal_places_ > 0) {
    DEBUG("decimal_places_ " << decimal_places_);
    const double max_abs = std::max(std::abs(u1), std::abs(u2));
    const double tol = std::pow(10, -decimal_places_);
    if (std::abs(u1 - u2)/max_abs > tol && max_abs > tol) {
      pass = false;
    }
  } else {
    if (std::abs(u1 - u2) > tolerance_) {
      pass = false;
    }
  }
  return pass;
}

std::string CheckEnergy::err_msg_() const {
  std::stringstream ss;
  if (decimal_places_ > 0) {
    ss << "number of decimal places(" << decimal_places_ << "). ";
  } else {
    ss << "tolerance(" << tolerance_ << "). ";
  }
  return ss.str();
}

void CheckEnergy::update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  check_->update(*criteria, *system, *trial_factory);
  DEBUG("computing unoptimized energy for check");

  for (int config = 0; config < system->num_configurations(); ++config) {
    DEBUG("config " << config);
    const double energy = system->unoptimized_energy(config);
    // HWH configuration_index_
    const double current_energy = criteria->current_energy(config);
    DEBUG("energy:" << energy << " "
       << "current_energy: " << current_energy << " "
       << "diff: " << energy - current_energy
    );
    // HWH configuration_index_
    accumulator_->accumulate(energy - current_energy);

    // loop over each profile and perform the energy check
    const std::vector<double>& energy_profile = system->unoptimized(config).stored_energy_profile();
    // HWH configuration_index_
    const std::vector<double>& current_energy_profile = criteria->current_energy_profile(config);
    DEBUG("energy_profile " << feasst_str(energy_profile));
    DEBUG("current_energy_profile " << feasst_str(current_energy_profile));
    for (int i = 0; i < static_cast<int>(energy_profile.size()); ++i) {
      ASSERT(is_within_tolerance_(energy_profile[i], current_energy_profile[i]),
        MAX_PRECISION <<
        "Energy check failure. There is a problem with the potentials. " <<
        "The unoptimized energy of potential " << i << " in configuration " <<
        config << " was computed as " <<
        energy_profile[i] << " but the running energy from criteria " <<
        "(the accumulation of a change in energy over a series of steps) is " <<
        current_energy_profile[i] <<
        ". The difference(" << std::abs(energy_profile[i] - current_energy_profile[i]) << ") is " <<
        "greater than the " << err_msg_());
    }
    criteria->set_current_energy_profile(energy_profile, config);

    // perform same energy check for entire system
    ASSERT(is_within_tolerance_(energy, current_energy),
      MAX_PRECISION <<
      "Energy check failure. There is a problem with the potentials. " <<
      "The unoptimized energy of configuration " << config <<
      " was computed as " << energy << " but the running energy from criteria "
      << "(the accumulation of a change in energy over a series of steps) is "
      << current_energy <<
      ". The difference(" << std::abs(energy - current_energy) << ") is " <<
      "greater than the " << err_msg_() << " "
      << system->unoptimized().str());
    criteria->set_current_energy(energy, config);

    // loop over all queryable maps and check those as well.
  }
}

void CheckEnergy::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(716, ostr);
  feasst_serialize(tolerance_, ostr);
  feasst_serialize(decimal_places_, ostr);
  feasst_serialize_fstdr(check_, ostr);
}

CheckEnergy::CheckEnergy(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 715 && version <= 716, "version mismatch: " << version);
  feasst_deserialize(&tolerance_, istr);
  if (version >= 716) {
    feasst_deserialize(&decimal_places_, istr);
  }
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
