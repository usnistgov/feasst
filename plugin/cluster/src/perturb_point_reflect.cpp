#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "cluster/include/perturb_point_reflect.h"

namespace feasst {

PerturbPointReflect::PerturbPointReflect(argtype args) : PerturbMove(&args) {
  class_name_ = "PerturbPointReflect";
  check_all_used(args);
}

class MapPerturbPointReflect {
 public:
  MapPerturbPointReflect() {
    auto obj = MakePerturbPointReflect();
    obj->deserialize_map()["PerturbPointReflect"] = obj;
  }
};

static MapPerturbPointReflect mapper_ = MapPerturbPointReflect();

void PerturbPointReflect::precompute(TrialSelect * select, System * system) {
  set_tunable_min_and_max(2*NEAR_ZERO,
    0.5*system->configuration().domain().max_side_length());
}

void PerturbPointReflect::update_selection(const Position& reflect,
    TrialSelect * select) {
  Select * reflected = select->get_mobile();
  for (int select_index = 0;
       select_index < reflected->num_particles();
       ++select_index) {
    for (int site = 0;
         site < static_cast<int>(reflected->site_indices(select_index).size());
         ++site) {
      reflected->get_site_position(select_index, site)->reflect(reflect);
    }
  }
}

void PerturbPointReflect::move(
    const Position& reflect,
    System * system,
    TrialSelect * select) {
  update_selection(reflect, select);
  system->get_configuration()->update_positions(select->mobile());
}

void PerturbPointReflect::move(System * system,
                            TrialSelect * select,
                            Random * random) {
  random->position_in_cube(
    system->dimension(),
    tunable().value(),
    &reflect_
  );
  ASSERT(select->mobile().num_particles() == 1, "assumes 1 particle");
  reflect_.add(select->mobile().site_positions()[0][0]);
  DEBUG("max move " << tunable().value());
  ASSERT(tunable().value() > NEAR_ZERO, "tunable(" << tunable().value()
    << ") is too small");
  move(reflect_, system, select);
}

std::shared_ptr<Perturb> PerturbPointReflect::create(std::istream& istr) const {
  return std::make_shared<PerturbPointReflect>(istr);
}

PerturbPointReflect::PerturbPointReflect(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbPointReflect", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(564 == version, "mismatch version: " << version);
}

void PerturbPointReflect::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(564, ostr);
}

}  // namespace feasst
