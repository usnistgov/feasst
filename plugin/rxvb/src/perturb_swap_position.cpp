#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "math/include/matrix.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "rxvb/include/perturb_swap_position.h"

namespace feasst {

PerturbSwapPosition::PerturbSwapPosition(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbSwapPosition";
  disable_tunable_();
}
PerturbSwapPosition::PerturbSwapPosition(argtype args) : PerturbSwapPosition(&args) {
  feasst_check_all_used(args);
}
PerturbSwapPosition::~PerturbSwapPosition() {}

FEASST_MAPPER(PerturbSwapPosition,);

std::shared_ptr<Perturb> PerturbSwapPosition::create(std::istream& istr) const {
  return std::make_shared<PerturbSwapPosition>(istr);
}

void PerturbSwapPosition::move(const bool is_position_held,
                           System * system,
                           TrialSelect * select,
                           Random * random,
                           Acceptance * acceptance) {
  if (is_position_held) return;
  Select * mobile = select->get_mobile();
  DEBUG("mobile " << mobile->str());
  ASSERT(mobile->num_particles() == 2, "Error");
  const int dimen = select->configuration(*system).dimension();
  if (!r0_) {
    //r0_ = std::make_unique<Position>(argtype({{"dimension", str(dimen)}}));
    //r1_ = std::make_unique<Position>(argtype({{"dimension", str(dimen)}}));
    //axis_ = std::make_unique<Position>(argtype({{"dimension", str(dimen)}}));
    //origin_ = std::make_unique<Position>(argtype({{"dimension", str(dimen)}}));
    quat_ = std::make_unique<Position>(dimen + 1);
    r0_ = std::make_unique<Position>(dimen);
    r1_ = std::make_unique<Position>(dimen);
    axis_ = std::make_unique<Position>(dimen);
    origin_ = std::make_unique<Position>(dimen);
    rot_mat_ = std::make_unique<RotationMatrix>();
  }
  *r0_ = mobile->site_positions()[0][0];
  *r1_ = mobile->site_positions()[1][0];

  // Translate particles to the origin
  r0_->multiply(-1);
  r1_->multiply(-1);
  for (int site = 0; site < mobile->num_sites(0); ++site) {
    mobile->add_to_site_position(0, site, *r0_);
  }
  for (int site = 0; site < mobile->num_sites(1); ++site) {
    mobile->add_to_site_position(1, site, *r1_);
  }

  // Rotate
  // Update PerturbRotate to work with only a portion of Select??
  // Otherwise need to worry about Eulers?
  ASSERT(select->is_isotropic(system), "Not implemented");
  random->rotation(dimen, quat_.get(), rot_mat_.get());
  for (int site = 0; site < mobile->num_sites(0); ++site) {
    rot_mat_->rotate(*origin_, mobile->get_site_position(0, site), axis_.get());
  }
  random->rotation(dimen, quat_.get(), rot_mat_.get());
  for (int site = 0; site < mobile->num_sites(1); ++site) {
    rot_mat_->rotate(*origin_, mobile->get_site_position(1, site), axis_.get());
  }

  // Translate particles to their new positions
  r0_->multiply(-1);
  r1_->multiply(-1);
  for (int site = 0; site < mobile->num_sites(0); ++site) {
    mobile->add_to_site_position(0, site, *r1_);
  }
  for (int site = 0; site < mobile->num_sites(1); ++site) {
    mobile->add_to_site_position(1, site, *r0_);
  }
  DEBUG("mobile " << mobile->str());
  select->get_configuration(system)->update_positions(select->mobile());
}

PerturbSwapPosition::PerturbSwapPosition(std::istream& istr)
  : PerturbMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(6036 == version, "mismatch version: " << version);
}

void PerturbSwapPosition::serialize_perturb_position_swap_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(6036, ostr);
}

void PerturbSwapPosition::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_position_swap_(ostr);
}

}  // namespace feasst
