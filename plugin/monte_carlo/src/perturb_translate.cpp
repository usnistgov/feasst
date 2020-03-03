
#include "monte_carlo/include/perturb_translate.h"

namespace feasst {

PerturbTranslate::PerturbTranslate(const argtype& args) : PerturbMove(args) {
  class_name_ = "PerturbTranslate";
}

class MapPerturbTranslate {
 public:
  MapPerturbTranslate() {
    auto obj = MakePerturbTranslate();
    obj->deserialize_map()["PerturbTranslate"] = obj;
  }
};

static MapPerturbTranslate mapper_ = MapPerturbTranslate();

void PerturbTranslate::precompute(TrialSelect * select, System * system) {
  set_tunable_min_and_max(2*NEAR_ZERO,
    0.5*system->configuration().domain()->max_side_length());
}

void PerturbTranslate::update_selection(const Position& trajectory,
    TrialSelect * select) {
  SelectList * displaced = select->get_mobile();
  for (int select_index = 0;
       select_index < displaced->num_particles();
       ++select_index) {
    displaced->add_to_particle_position(select_index, trajectory);
    for (int site = 0;
         site < static_cast<int>(displaced->site_indices(select_index).size());
         ++site) {
      displaced->add_to_site_position(select_index, site, trajectory);
    }
  }
}

void PerturbTranslate::move(
    const Position& trajectory,
    System * system,
    TrialSelect * select) {
  update_selection(trajectory, select);
  system->get_configuration()->update_positions(select->mobile());
}

void PerturbTranslate::move(System * system,
                            TrialSelect * select,
                            Random * random) {
  random->position_in_cube(
    system->dimension(),
    tunable().value(),
    &trajectory_
  );
  DEBUG("max move " << tunable().value());
  ASSERT(tunable().value() > NEAR_ZERO, "tunable(" << tunable().value()
    << ") is too small");
  move(trajectory_, system, select);
}

std::shared_ptr<Perturb> PerturbTranslate::create(std::istream& istr) const {
  return std::make_shared<PerturbTranslate>(istr);
}

PerturbTranslate::PerturbTranslate(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbTranslate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(564 == version, "mismatch version: " << version);
}

void PerturbTranslate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(564, ostr);
}

}  // namespace feasst
