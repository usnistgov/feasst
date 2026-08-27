
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "rxvb/include/select_two_particles.h"

namespace feasst {

FEASST_MAPPER(SelectTwoParticles,);

SelectTwoParticles::SelectTwoParticles(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectTwoParticles";
  // HWH update this to comma-separated in v0.25.13
  pt_names_ = split(str("particle_types", args, ""), ',');
  const int num = static_cast<int>(pt_names_.size());
  ASSERT(num == 0 || num == 2,
    "Unrecognized format for particle_types argument");
  for (int ipt1 = 0; ipt1 < static_cast<int>(pt_names_.size())-1; ++ipt1) {
    for (int ipt2 = ipt1+1; ipt2 < static_cast<int>(pt_names_.size()); ++ipt2) {
      ASSERT(pt_names_[ipt1] != pt_names_[ipt2], "SelectTwoParticles requires that" <<
        " particle_types lists different particles. But particle type \"" <<
        pt_names_[ipt1] << "\" appears more than once.");
    }
  }
}
SelectTwoParticles::SelectTwoParticles(argtype args) : SelectTwoParticles(&args) {
  feasst_check_all_used(args);
}
SelectTwoParticles::~SelectTwoParticles() {}

std::shared_ptr<TrialSelect> SelectTwoParticles::create(std::istream& istr) const {
  return std::make_shared<SelectTwoParticles>(istr);
}

SelectTwoParticles::SelectTwoParticles(std::istream& istr) : TrialSelect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 2743 && version <= 2743, "mismatch version: " << version);
  feasst_deserialize(&pts_, istr);
  feasst_deserialize(&pt_names_, istr);
  feasst_deserialize(&group_indices_, istr);
}

void SelectTwoParticles::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_(ostr);
  feasst_serialize_version(2743, ostr);
  feasst_serialize(pts_, ostr);
  feasst_serialize(pt_names_, ostr);
  feasst_serialize(group_indices_, ostr);
}

void SelectTwoParticles::precompute(System * system) {
  TrialSelect::precompute(system);
  ASSERT(static_cast<int>(pt_names_.size()) == 2, "particle_types not provided.");
  group_indices_.resize(2);
  pts_.resize(2);
  for (int ipt = 0; ipt < static_cast<int>(group_indices_.size()); ++ipt) {
    const Configuration& config = configuration(*system);
    pts_[ipt] = config.particle_name_to_type(pt_names_[ipt]);
    group_indices_[ipt] =
      get_configuration(system)->particle_type_to_group_create(pts_[ipt]);
  }
  DEBUG("grp:" << feasst_str(group_indices_));
}

bool SelectTwoParticles::select(const Select& perturbed,
    System* system,
    Random * random,
    TrialSelect * previous_select) {
  const Configuration& config = configuration(*system);
  ASSERT(static_cast<int>(group_indices_.size()) == 2, "err");
  const int num0 = config.num_particles(group_indices_[0]);
  if (num0 <= 0) return false;
  const int p0index = random->uniform(0, num0 - 1);
  const Select& grp0 = config.group_select(group_indices_[0]);
  const int num1 = config.num_particles(group_indices_[1]);
  if (num1 <= 0) return false;
  const int p1index = random->uniform(0, num1 - 1);
  const Select& grp1 = config.group_select(group_indices_[1]);
  const bool fast = get_mobile()->replace_indices(
    {grp0.particle_index(p0index), grp1.particle_index(p1index)},
    {grp0.site_indices(p0index), grp1.site_indices(p1index)});
  if (!fast) get_mobile()->resize_positions();
  get_mobile()->load_positions(config.particles());
  remove_unphysical_sites(config);
  set_mobile_original(system);
  return true;
}

}  // namespace feasst
