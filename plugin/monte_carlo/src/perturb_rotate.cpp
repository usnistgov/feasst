#include <cmath>
#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "math/include/constants.h"
#include "math/include/random.h"

namespace feasst {

PerturbRotate::PerturbRotate(const argtype& args) : PerturbMove(args) {
  class_name_ = "PerturbRotate";
  set_tunable_min_and_max(2*NEAR_ZERO, 360.);
  Arguments args_(args);
  args_.dont_check();
  pivot_site_ = args_.key("pivot_site").dflt("-1").integer();
}

class MapPerturbRotate {
 public:
  MapPerturbRotate() {
    auto obj = MakePerturbRotate();
    obj->deserialize_map()["PerturbRotate"] = obj;
  }
};

static MapPerturbRotate mapper_ = MapPerturbRotate();

std::shared_ptr<Perturb> PerturbRotate::create(std::istream& istr) const {
  return std::make_shared<PerturbRotate>(istr);
}

void PerturbRotate::update_selection(const Position& pivot,
    const RotationMatrix& rotation,
    TrialSelect * select,
    const bool rotate_particle_position) {
  Select * rotated = select->get_mobile();
  for (int select_index = 0;
       select_index < rotated->num_particles();
       ++select_index) {
    // rotate site positions
    for (int site = 0;
         site < static_cast<int>(rotated->site_indices(select_index).size());
         ++site) {
      Position position = rotated->site_positions()[select_index][site];
      rotation.rotate(pivot, &position);
      rotated->set_site_position(select_index, site, position);
    }

    // rotate or recenter particle positions
    if (rotate_particle_position) {
      Position position = rotated->particle_positions()[select_index];
      rotation.rotate(pivot, &position);
      rotated->set_particle_position(select_index, position);
    }
  }
}

void PerturbRotate::move(const Position& pivot,
    const RotationMatrix& rotation,
    System * system,
    TrialSelect * select,
    const bool rotate_particle_position) {
  update_selection(pivot, rotation, select, rotate_particle_position);
  system->get_configuration()->update_positions(select->mobile());
}

void PerturbRotate::move(System * system,
    TrialSelect * select,
    Random * random) {
  ASSERT(select->mobile().num_sites() > 0, "selection error");
  const Position * pivot;
  if (pivot_site_ == -1) {
    pivot = const_cast<Position *>(&select->mobile().particle_positions()[0]);
  } else {
    pivot = const_cast<Position *>(
      &select->mobile().site_positions()[0][pivot_site_]);
  }
  move(system, select, random, *pivot, true);
}

void PerturbRotate::move(System * system,
    TrialSelect * select,
    Random * random,
    const Position& pivot,
    const bool rotate_particle_position) {
  if (is_rotation_not_needed_(select, pivot)) return;
  const double max_angle = tunable().value();
  ASSERT(std::abs(max_angle) > NEAR_ZERO, "max angle is too small");
  const Position& piv_sel = piv_sel_(pivot, select);
  random->rotation(piv_sel.dimension(), &axis_tmp_, &rot_mat_tmp_, max_angle),
  move(piv_sel,
    rot_mat_tmp_,
    system,
    select,
    rotate_particle_position
  );
}

PerturbRotate::PerturbRotate(std::istream& istr)
  : PerturbMove(istr) {
  // ASSERT(class_name_ == "PerturbRotate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(448 == version, "mismatch version: " << version);
  feasst_deserialize(&pivot_site_, istr);
}

void PerturbRotate::serialize_perturb_rotate_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(448, ostr);
  feasst_serialize(pivot_site_, ostr);
}

void PerturbRotate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_rotate_(ostr);
}

}  // namespace feasst
