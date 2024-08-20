#include <fstream>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/formula.h"
#include "math/include/golden_search.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/monte_carlo.h"
#include "patch/include/two_particle_contact.h"

namespace feasst {

TwoParticleContact::TwoParticleContact(argtype * args) {
  class_name_ = "TwoParticleContact";
  output_file_ = str("output_file", args, "");
  dimension_ = integer("dimension", args, 0);
  tolerance_ = dble("tolerance", args, 1e-6);
}
TwoParticleContact::TwoParticleContact(argtype args) : TwoParticleContact(&args) {
  feasst_check_all_used(args);
}

class MapTwoParticleContact {
 public:
  MapTwoParticleContact() {
    auto obj = MakeTwoParticleContact();
    obj->deserialize_map()["TwoParticleContact"] = obj;
  }
};

static MapTwoParticleContact mapper_TwoParticleContact = MapTwoParticleContact();

TwoParticleContact::TwoParticleContact(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3467 && version <= 3467, "mismatch version: " << version);
  feasst_deserialize(&output_file_, istr);
  feasst_deserialize(&dimension_, istr);
  feasst_deserialize(&tolerance_, istr);
}

void TwoParticleContact::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3467, ostr);
  feasst_serialize(output_file_, ostr);
  feasst_serialize(dimension_, ostr);
  feasst_serialize(tolerance_, ostr);
}

void TwoParticleContact::run(MonteCarlo * mc) {
  DEBUG("TwoParticleContact");
  ASSERT(mc, "mc required");
  System * system = mc->get_system();
  select_ = MakeTrialSelectParticle({{"particle_type", "1"}});
  select_->select_particle(1, system->configuration());
  select_->set_mobile_original(system);
  translate_ = MakePerturbTranslate();
  GoldenSearch minimize({{"tolerance", str(tolerance_)},
    {"lower", "0"},
    {"upper", str(system->configuration().domain().side_length(dimension_)/2)}});
  TwoParticleContactObjective objective(this, system);
  contact_distance_ = minimize.minimum(&objective) + 2.*tolerance_;
  if (!output_file_.empty()) {
    std::ofstream file(output_file_);
    ASSERT(file.good(), "err");
    file << contact_distance_;
    file.close();
  }
}

TwoParticleContactObjective::TwoParticleContactObjective(TwoParticleContact * two_particle_contact, System * system) {
  two_particle_contact_ = two_particle_contact;
  system_ = system;
}

double TwoParticleContactObjective::evaluate(const double distance) const {
  const double en = two_particle_contact_->energy(distance, system_);
  DEBUG("dist " << distance << " en " << en);
  if (en > 1) {
    return 1e15/(distance + 0.1);
  } else {
    return distance + 0.1;
  }
}

double TwoParticleContact::energy(const double distance, System * system) {
  update_xyz(distance, system);
  double en;
  //ASSERT(select_->mobile().num_particles() == 1, "err");
  en = system->perturbed_energy(select_->mobile());
  revert(system);
  DEBUG("position of first site of mobile after revert " << system->configuration().particle(1).site(0).position().str());
  return en;
}

void TwoParticleContact::update_xyz(const double distance, System * system) {
  com1_.set_to_origin(3);
  com1_.set_coord(dimension_, distance);
  DEBUG("com1 " << com1_.str());
  select_->select_particle(1, system->configuration());
  translate_->set_revert_possible(true, select_.get());
  translate_->move(com1_, system, select_.get());
  DEBUG("first atom " << system->configuration().particle(1).site(0).position().str());
  DEBUG("first atom fixed " << system->configuration().particle(0).site(0).position().str());
}

void TwoParticleContact::revert(System * system) {
  DEBUG("position of first site of mobile before revert " << system->configuration().particle(1).site(0).position().str());
  system->revert(select_->mobile());
  translate_->revert(system);
  DEBUG("position of first site of mobile after revert " << system->configuration().particle(1).site(0).position().str());
}

}  // namespace feasst
