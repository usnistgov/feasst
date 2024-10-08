#include <cmath>
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "system/include/thermo_params.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "steppers/include/heat_capacity.h"

namespace feasst {

FEASST_MAPPER(HeatCapacity,);

HeatCapacity::HeatCapacity(argtype * args) : Analyze(args) {
  energy_ = *MakeAccumulator({{"num_moments", "3"}});
}
HeatCapacity::HeatCapacity(argtype args) : HeatCapacity(&args) {
  feasst_check_all_used(args);
}

void HeatCapacity::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          output_file(*criteria));
}

std::string HeatCapacity::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << "heat_capacity_per_kB" << std::endl;
  return ss.str();
}

void HeatCapacity::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const double en = criteria.current_energy(configuration_index());
  DEBUG("en: " << en);
  DEBUG("state: " << state());
  energy_.accumulate(en);
}

std::string HeatCapacity::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(criteria, system, trial_factory);
  }
  const double beta = system.thermo_params().beta();
  const double u_sq_av = energy_.moment(2)/energy_.moment(0);
  const double u_av_sq = std::pow(energy_.average(), 2);
  ss << beta*beta*(u_sq_av - u_av_sq) << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void HeatCapacity::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2347, ostr);
  feasst_serialize_fstobj(energy_, ostr);
}

HeatCapacity::HeatCapacity(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2347, "mismatch version:" << version);
  feasst_deserialize_fstobj(&energy_, istr);
}

HeatCapacity::HeatCapacity(const Analyze& energy) {
  std::stringstream ss;
  energy.serialize(ss);
  *this = HeatCapacity(ss);
}

}  // namespace feasst
