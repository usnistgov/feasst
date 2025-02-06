#include <cmath>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "monte_carlo/include/tunable.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/perturb_volume.h"

namespace feasst {

PerturbVolume::PerturbVolume(argtype * args) : Perturb(args) {
  class_name_ = "PerturbVolume";
  uniform_volume_ = boolean("uniform_volume", args, false);
  constrain_volume_change_ = boolean("constrain_volume_change", args, false);
  if (constrain_volume_change_) {
    disable_tunable_();
  }
}
PerturbVolume::PerturbVolume(argtype args) : PerturbVolume(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(PerturbVolume,);

std::shared_ptr<Perturb> PerturbVolume::create(std::istream& istr) const {
  return std::make_shared<PerturbVolume>(istr);
}

void PerturbVolume::precompute(TrialSelect * select, System * system) {
  const int config = select->configuration_index();
  args_.insert({"configuration", str(config)});
}

void PerturbVolume::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held,
    Acceptance * acceptance) {
  ASSERT(!is_position_held, "not implemeted");
  DEBUG("config " << select->configuration_index());
  const double volume = select->configuration(*system).domain().volume();
  DEBUG("volume " << MAX_PRECISION << volume);
  if (constrain_volume_change_) {
    volume_change_ = -system->delta_volume_previous();
  } else  if (uniform_volume_) {
    volume_change_ = random->uniform_real(-tunable().value(),
                                           tunable().value());
  } else {
    // lnvn = lnv0 + dlnv
    //vn = exp(lnvo + dlnv)
    //dv = exp(lnvo + dlnv) - vo
    const double dlnv = random->uniform_real(-tunable().value(),
                                              tunable().value());
    volume_change_ = std::exp(std::log(volume) + dlnv) - volume;
  }
  DEBUG("volume_change_ " << volume_change_);
  if (volume + volume_change_ > 0) {
    change_volume(volume_change_, system, select->mobile());
  } else {
    volume_change_ = 0.;
  }
  DEBUG("volume_change_ " << volume_change_);
  set_revert_possible(true, select);
  set_finalize_possible(true, select);
  select->set_trial_state(4);
}

void PerturbVolume::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
    DEBUG(revert_select()->mobile().str());
    DEBUG("nump " << system->configuration().num_particles());
//    system->revert(revert_select()->mobile());
    DEBUG("volume_change_ " << volume_change_);
    change_volume(-volume_change_, system, revert_select()->mobile());
  }
}

void PerturbVolume::finalize(System * system) {
  DEBUG("finalize_possible " << finalize_possible());
  if (finalize_possible()) {
//    system->finalize(finalize_select()->mobile());
  }
}

PerturbVolume::PerturbVolume(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbVolume", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 1634 && version <= 1635, "mismatch version: " << version);
  if (version >= 1635) {
    feasst_deserialize(&uniform_volume_, istr);
    feasst_deserialize(&constrain_volume_change_, istr);
    feasst_deserialize(&args_, istr);
  }
}

void PerturbVolume::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(1635, ostr);
  feasst_serialize(uniform_volume_, ostr);
  feasst_serialize(constrain_volume_change_, ostr);
  feasst_serialize(args_, ostr);
}

void PerturbVolume::change_volume(const double delta_volume,
    System * system,
    const Select& select) {
  system->change_volume(delta_volume, args_);
}

}  // namespace feasst
