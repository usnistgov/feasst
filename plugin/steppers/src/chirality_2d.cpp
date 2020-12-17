#include <cmath>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "steppers/include/chirality_2d.h"

namespace feasst {

class MapChirality2D {
 public:
  MapChirality2D() {
    Chirality2D().deserialize_map()["Chirality2D"] = MakeChirality2D();
  }
};

static MapChirality2D mapper_energy_check_ = MapChirality2D();

Chirality2D::Chirality2D(const argtype &args)
  : Analyze(args) {
  group_ = args_.key("group").dflt("0").integer();
  bond1_ = args_.key("bond1").dflt("0").integer();
  bond2_ = args_.key("bond2").dflt("1").integer();
  sign_error_ = args_.key("sign_error").dflt("0").integer();
}

std::string Chirality2D::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator_.status_header() << std::endl;
  return ss.str();
}

void Chirality2D::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  ASSERT(system->configuration().dimension() == 2,
    "dim: " << system->configuration().dimension() << " != 2");
  printer(header(*criteria, *system, *trial_factory),
    file_name(*criteria));
}

void Chirality2D::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
    const Configuration& config = system.configuration();
  int num_positive = 0;
  for (int part_index : config.group_selects()[0].particle_indices()) {
    const Particle& part = config.select_particle(part_index);
    const Bond& bond1 = config.particle_type(part.type()).bond(bond1_);
    const Bond& bond2 = config.particle_type(part.type()).bond(bond2_);
    Position pos1(3), pos2(3);
    for (int dim = 0; dim < config.dimension(); ++dim) {
      pos1.set_coord(dim, part.site(bond1.site(1)).position().coord(dim)
                        - part.site(bond1.site(0)).position().coord(dim));
      pos2.set_coord(dim, part.site(bond2.site(1)).position().coord(dim)
                        - part.site(bond2.site(0)).position().coord(dim));
    }
    const double cross_z = pos1.cross_product(pos2).coord(2);
    if (cross_z) ++num_positive;
    if (sign_error_ != 0) {
      if (sign_error_ > 0 && cross_z > 0) FATAL("positive chirality");
      if (sign_error_ < 0 && cross_z < 0) FATAL("negative chirality");
    }
  }
  accumulator_.accumulate(num_positive);
}

std::string Chirality2D::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << accumulator_.status() << std::endl;
  return ss.str();
}

void Chirality2D::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(732, ostr);
  feasst_serialize(group_, ostr);
  feasst_serialize(bond1_, ostr);
  feasst_serialize(bond2_, ostr);
  feasst_serialize(sign_error_, ostr);
}

Chirality2D::Chirality2D(std::istream& istr)
  : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(732 == version, "version mismatch:" << version);
  feasst_deserialize(&group_, istr);
  feasst_deserialize(&bond1_, istr);
  feasst_deserialize(&bond2_, istr);
  feasst_deserialize(&sign_error_, istr);
}

}  // namespace feasst
