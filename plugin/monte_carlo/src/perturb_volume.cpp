#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/perturb_volume.h"

namespace feasst {

PerturbVolume::PerturbVolume(const argtype& args) : Perturb(args) {
  class_name_ = "PerturbVolume";
}

class MapPerturbVolume {
 public:
  MapPerturbVolume() {
    auto obj = MakePerturbVolume();
    obj->deserialize_map()["PerturbVolume"] = obj;
  }
};

static MapPerturbVolume mapper_ = MapPerturbVolume();

std::shared_ptr<Perturb> PerturbVolume::create(std::istream& istr) const {
  return std::make_shared<PerturbVolume>(istr);
}

void PerturbVolume::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held) {
  ASSERT(!is_position_held, "not implemented");
  // lnvn = lnv0 + dlnv
  //vn = exp(lnvo + dlnv)
  //dv = exp(lnvo + dlnv) - vo
  const double dlnv = random->uniform_real(-tunable().value(),
                                            tunable().value());
  const double volume = system->configuration().domain().volume();
  volume_change_ = std::exp(std::log(volume) + dlnv) - volume;
  //volume_change_ = dlnv;  // temporary test v' = v + dv
  if (volume + volume_change_ > 0) {
    change_volume(volume_change_, system, select->mobile());
  } else {
    volume_change_ = 0.;
  }
  set_revert_possible(true, select);
  set_finalize_possible(true, select);
  select->set_trial_state(1);
}

void PerturbVolume::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
    DEBUG(revert_select()->mobile().str());
    DEBUG("nump " << system->configuration().num_particles());
    system->revert(revert_select()->mobile());
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
  ASSERT(1634 == version, "mismatch version: " << version);
}

void PerturbVolume::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(1634, ostr);
}

}  // namespace feasst
