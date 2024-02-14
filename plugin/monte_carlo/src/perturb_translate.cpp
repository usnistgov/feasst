#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "monte_carlo/include/perturb_translate.h"

namespace feasst {

PerturbTranslate::PerturbTranslate(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbTranslate";
  dimension_ = integer("dimension", args, -1);
  if (dimension_ != -1) {
    disable_tunable_();
  }
}
PerturbTranslate::PerturbTranslate(argtype args) : PerturbTranslate(&args) {
  FEASST_CHECK_ALL_USED(args);
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
  const Configuration& config = select->configuration(*system);
  set_tunable_min_and_max(2*NEAR_ZERO,
    0.5*config.domain().max_side_length());  // +/- L/2
}

void PerturbTranslate::update_selection(const Position& trajectory,
    TrialSelect * select) {
  DEBUG("traj " << trajectory.str());
  Select * displaced = select->get_mobile();
  DEBUG("mob b4 " << displaced->site_positions()[0][0].str());
  if (anchor_set_) {
    new_pos_.set_to_origin(anchor_.size());
    new_pos_.add(anchor_);
    new_pos_.add(trajectory);
  }
  for (int select_index = 0;
       select_index < displaced->num_particles();
       ++select_index) {
    for (int site = 0;
         site < static_cast<int>(displaced->site_indices(select_index).size());
         ++site) {
      if (anchor_set_) {
        displaced->set_site_position(select_index, site, new_pos_);
      } else {
        displaced->add_to_site_position(select_index, site, trajectory);
      }
    }
  }
  DEBUG("mob af " << displaced->site_positions()[0][0].str());
}

void PerturbTranslate::move(
    const Position& trajectory,
    System * system,
    TrialSelect * select) {
  DEBUG("trajectory " << trajectory.str());
  update_selection(trajectory, select);
  select->get_configuration(system)->update_positions(select->mobile());
}

void PerturbTranslate::move(const bool is_position_held,
                            System * system,
                            TrialSelect * select,
                            Random * random) {
  if (is_position_held) return;
  const int dim = system->dimension();
  if (dimension_ == -1) {
    random->position_in_cube(
      dim,
      2.*tunable().value(), // 2 factor accounts for +/- delta
      &trajectory_
    );
  } else {
    ASSERT(dimension_ < dim, "dimension_: " << dimension_ << " < " << dim);
    ASSERT(!tunable().is_enabled(),
      "PeturbTranslate assumes tunable is disabled.");
    trajectory_.set_to_origin(dim);
    int sign = 1.;
    if (random->coin_flip()) {
      sign = -1;
    }
    trajectory_.set_coord(dimension_, sign*tunable().value());
  }
  DEBUG("max move " << tunable().value());
  move(trajectory_, system, select);
}

void PerturbTranslate::begin_stage(const TrialSelect& select) {
  if (select.mobile().num_sites() == 1) {
    anchor_ = select.mobile().site_positions()[0][0];
    anchor_set_ = true;
  } else {
    anchor_set_ = false;
  }
}

void PerturbTranslate::mid_stage(const TrialSelect& select,
    const System& system) {
  if (select.mobile().num_sites() == 1) {
    // Fix this , mobile isn't in the chosen position!
    //anchor_ = select.mobile().site_positions()[0][0];
    anchor_ = select.configuration(system).select_particle(select.mobile().particle_index(0)).site(select.mobile().site_index(0, 0)).position();
    anchor_set_ = true;
  }
}

std::shared_ptr<Perturb> PerturbTranslate::create(std::istream& istr) const {
  return std::make_shared<PerturbTranslate>(istr);
}

PerturbTranslate::PerturbTranslate(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbTranslate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 564 && version <= 566, "mismatch version: " << version);
  if (version <= 565) {
    // In previous version, erroneously serialized anchor_set_
    bool temporary;
    feasst_deserialize(&temporary, istr);
  }
  if (version >= 565) {
    feasst_deserialize(&dimension_, istr);
  }
}

void PerturbTranslate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(566, ostr);
  feasst_serialize(dimension_, ostr);
}

}  // namespace feasst
