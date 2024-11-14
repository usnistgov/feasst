#include <cmath>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/formula.h"
#include "math/include/accumulator.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/criteria.h"
#include "steppers/include/pressure_from_test_volume.h"

namespace feasst {

FEASST_MAPPER(PressureFromTestVolume,);

PressureFromTestVolume::~PressureFromTestVolume() {}

void PressureFromTestVolume::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(8947, ostr);
  feasst_serialize(delta_volume_, ostr);
}

PressureFromTestVolume::PressureFromTestVolume(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(8947 == version, version);
  feasst_deserialize(&delta_volume_, istr);
}

PressureFromTestVolume::PressureFromTestVolume(argtype * args) : Modify(args) {
  delta_volume_ = dble("delta_volume", args, 1e-4);
}
PressureFromTestVolume::PressureFromTestVolume(argtype args) : PressureFromTestVolume(&args) {
  feasst_check_all_used(args);
}

void PressureFromTestVolume::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
}

void PressureFromTestVolume::update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  DEBUG("updating");
  const int config_id = configuration_index();
  const double en_old = criteria->current_energy(config_id);
  argtype args;
  args.insert({"configuration", str(config_id)});
  system->change_volume(delta_volume_, args);
  const double en_new = system->energy(config_id);
  system->change_volume(-delta_volume_, args);
  const double volume = configuration(*system).domain().volume();
  const double num_particles = configuration(*system).num_particles();
  const double beta = system->thermo_params().beta();
  const double ens_av = std::pow((volume + delta_volume_)/volume, num_particles)*std::exp(-beta*(en_new - en_old))-1.;
  DEBUG("ens_av " << MAX_PRECISION << ens_av);
  get_accumulator()->accumulate(ens_av);
}

std::string PressureFromTestVolume::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << "average,block_stdev,";
  for (int index = 0;
       index < static_cast<int>(accumulator().largest_blocks().size());
       ++index) {
    ss << "block" << index << ",";
  }
  ss << std::endl;
  return ss.str();
}

class LogEnsAv : public Formula {
 public:
  double evaluate(const double y) const override { return std::log(y+1.); }
};


std::string PressureFromTestVolume::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  std::stringstream ss;
  ss << header(*criteria, *system, *trial_factory);
  const double beta = system->thermo_params().beta();
  LogEnsAv les;
  ss << std::log(accumulator().average()+1.)/beta/delta_volume_ << ","
     << accumulator().block_stdev(les)/beta/delta_volume_ << ",";
  for (double val : accumulator().largest_blocks()) {
    ss << std::log(val+1.)/beta/delta_volume_ << ",";
  }
  return ss.str();
}

}  // namespace feasst
