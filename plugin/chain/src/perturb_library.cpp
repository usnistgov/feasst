#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/trial_select.h"
#include "chain/include/perturb_library.h"

namespace feasst {

FEASST_MAPPER(PerturbLibrary, argtype({{"library_xyz", "empty"}}));

PerturbLibrary::PerturbLibrary(argtype args) : PerturbLibrary(&args) {
  feasst_check_all_used(args);
}
PerturbLibrary::PerturbLibrary(argtype * args) : PerturbRotate(args) {
  class_name_ = "PerturbLibrary";
  library_xyz_file_name_ = str("library_xyz", args);
}

std::shared_ptr<Perturb> PerturbLibrary::create(std::istream& istr) const {
  return std::make_shared<PerturbLibrary>(istr);
}

PerturbLibrary::PerturbLibrary(std::istream& istr)
  : PerturbRotate(istr) {
  ASSERT(class_name_ == "PerturbLibrary", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(6079 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&xyz_, istr);
}

void PerturbLibrary::serialize_perturb_library_(std::ostream& ostr) const {
  serialize_perturb_rotate_(ostr);
  feasst_serialize_version(6079, ostr);
  feasst_serialize_fstobj(xyz_, ostr);
}

void PerturbLibrary::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_library_(ostr);
}

void PerturbLibrary::precompute(TrialSelect * select, System * system) {
  FileXYZ file_xyz;
  std::ifstream xyz(library_xyz_file_name_);
  ASSERT(xyz, library_xyz_file_name_ << " is empty");
  const Configuration& config = select->configuration(*system);
  auto one_particle_config = MakeConfiguration({
    {"particle_type0", config.type_to_file_name(select->particle_type())}});
  const int num_sites = one_particle_config->particle_type(0).num_sites();
  std::vector<Position> one_xyz(num_sites);
  while (xyz.peek() != EOF) {
    file_xyz.load_frame(xyz, one_particle_config.get());
    const Particle& part = one_particle_config->particle(0);
    for (int site = 0; site < num_sites; ++site) {
      one_xyz[site] = part.site(site).position();
    }
    xyz_.push_back(one_xyz);
  }
}

void PerturbLibrary::move(const bool is_position_held,
                        System * system,
                        TrialSelect * select,
                        Random * random,
                        Acceptance * acceptance) {
  FATAL("not implemented");
//  if (is_position_held) return;
  //const Configuration& config = system->configuration();
  //ASSERT(site == -1, "PerturbLibrary requires whole particle");
  //const int pivot = ;

//  const Position& pivot = select->anchor_position(0, 0, *system);
//  DEBUG("piv " << pivot.str());
//  PerturbRotate::move(system, select, random, pivot);
//  DEBUG(select->mobile().site_positions()[0][0].str());
}

}  // namespace feasst
