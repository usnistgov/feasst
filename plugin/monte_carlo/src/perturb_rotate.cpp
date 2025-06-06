#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/tunable.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

PerturbRotate::PerturbRotate(argtype args) : PerturbRotate(&args) {
  feasst_check_all_used(args);
}
PerturbRotate::PerturbRotate(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbRotate";
  set_tunable_min_and_max(2*NEAR_ZERO, 180.*(1 + NEAR_ZERO));
//  pivot_site_ = integer("pivot_site", args, 0);
}

FEASST_MAPPER(PerturbRotate,);

std::shared_ptr<Perturb> PerturbRotate::create(std::istream& istr) const {
  return std::make_shared<PerturbRotate>(istr);
}

void PerturbRotate::update_selection(const Position& pivot,
    const RotationMatrix& rotation,
    TrialSelect * select) {
  Select * rotated = select->get_mobile();
  DEBUG("rotation matrix : " << rotation.str());
  DEBUG("rotated " << rotated->str());
  for (int select_index = 0;
       select_index < rotated->num_particles();
       ++select_index) {
    // rotate site positions
    for (int site = 0;
         site < static_cast<int>(rotated->site_indices(select_index).size());
         ++site) {
      Position position = rotated->site_positions()[select_index][site];
      //rotation.rotate(pivot, &position);
      DEBUG("site " << site);
      DEBUG("old position " << position.str());
      rotation.rotate(pivot, &position, &vec2_);
      DEBUG("new position " << position.str());
      rotated->set_site_position(select_index, site, position);
    }
  }
}

void PerturbRotate::update_eulers(const RotationMatrix& rotation,
    TrialSelect * select,
    const System * system) {
  if (select->is_isotropic(system)) {
    return;
  }
  DEBUG("rotation " << rotation.str());
  const Configuration& config = select->configuration(*system);
  Select * rotated = select->get_mobile();
  for (int select_index = 0;
       select_index < rotated->num_particles();
       ++select_index) {
    const int part_index = rotated->particle_index(select_index);
    for (int select_site = 0;
         select_site < static_cast<int>(rotated->site_indices(select_index).size());
         ++select_site) {
      const int site_index = rotated->site_index(select_index, select_site);
//      DEBUG("part_index " << part_index);
      const Particle& part = config.select_particle(part_index);
//      DEBUG("site_index " << site_index);
      const Site& site = part.site(site_index);
      const int type = site.type();
//      DEBUG("type " << type);
      if (config.unique_type(part.type(), type).is_anisotropic()) {
        site.euler().compute_rotation_matrix(&rot_mat3_);
        DEBUG("rotation " << rotation.str());
        rotation.multiply(rot_mat3_, &rot_mat2_, &axis_tmp_, &vec1_);
        DEBUG("rot_mat3_ " << rot_mat3_.str());
        DEBUG("rot_mat2_ " << rot_mat2_.str());
        euler_.set(rot_mat2_);
        DEBUG("select_index " << select_index);
        DEBUG("select_site " << select_site);
        rotated->set_euler(select_index, select_site, euler_);
        DEBUG("new Euler(" << part_index << "," << site_index << ") " << euler_.str());
      }
    }
  }
}

void PerturbRotate::move(const Position& pivot,
    const RotationMatrix& rotation,
    System * system,
    TrialSelect * select) {
  update_selection(pivot, rotation, select);
  update_eulers(rotation, select, system);
  select->get_configuration(system)->update_positions(select->mobile());
}

void PerturbRotate::move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) {
  if (is_position_held) return;
  ASSERT(select->mobile().num_sites() > 0, "selection error");
  ASSERT(select->mobile().site_positions().size() > 0, "requires coordinates");
  const Position& pivot = select->mobile().site_positions()[0][0];
//  const Position& pivot = select->mobile().site_positions()[0][pivot_site_];
  move(system, select, random, pivot);
}

void PerturbRotate::move(System * system,
    TrialSelect * select,
    Random * random,
    const Position& pivot) {
  if (is_rotation_not_needed_(select, pivot)) {
    if (select->is_isotropic(system)) {
      return;
    }
  }
  const double tune = tunable().value();
  ASSERT(std::abs(tune) > NEAR_ZERO, "max angle is too small");
  const Position& piv_sel = piv_sel_(pivot, select);
  random->rotation(piv_sel.dimension(), &axis_tmp_, &rot_mat1_, tune);
  move(piv_sel, rot_mat1_, system, select);
}

PerturbRotate::PerturbRotate(std::istream& istr)
  : PerturbMove(istr) {
  // ASSERT(class_name_ == "PerturbRotate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 448 && version <= 449, "mismatch version: " << version);
  if (version == 448) {
    int temporary_value;
    feasst_deserialize(&temporary_value, istr);
  }
}

void PerturbRotate::serialize_perturb_rotate_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(449, ostr);
//  feasst_serialize(pivot_site_, ostr);
}

void PerturbRotate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_rotate_(ostr);
}

const Position& PerturbRotate::piv_sel_(const Position& pivot,
                                        const TrialSelect * select) {
  if (pivot.dimension() == 0) {
    return select->mobile().site_positions()[0][0];
  }
  return pivot;
}

bool PerturbRotate::is_rotation_not_needed_(const TrialSelect * select,
    const Position& pivot) {
  const Select& rotated = select->mobile();
  if (rotated.num_particles() == 1) {
    if (static_cast<int>(rotated.site_indices()[0].size()) == 1) {
      if (rotated.site_positions()[0][0].is_equal(pivot)) {
        return true;
      }
    }
  }
  return false;
}

}  // namespace feasst
