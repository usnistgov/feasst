#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/chirality_2d.h"

namespace feasst {

FEASST_MAPPER(Chirality2D,);

Chirality2D::Chirality2D(argtype * args) : Analyze(args) {
  group_ = integer("group", args, 0);
  bond1_ = integer("bond1", args, 0);
  bond2_ = integer("bond2", args, 1);
  sign_error_ = integer("sign_error", args, 0);
}
Chirality2D::Chirality2D(argtype args) : Chirality2D(&args) {
  feasst_check_all_used(args);
}

std::string Chirality2D::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator_->status_header() << std::endl;
  return ss.str();
}

void Chirality2D::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  System * system = mc->get_system();
  ASSERT(system->configuration().dimension() == 2,
    "dim: " << system->configuration().dimension() << " != 2");
  printer(header(*mc), output_file(mc->criteria()));
}

void Chirality2D::update(const MonteCarlo& mc) {
  const Configuration& config = mc.system().configuration();
  int num_positive = 0;
  for (int part_index : config.group_select(0).particle_indices()) {
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
  accumulator_->accumulate(num_positive);
}

std::string Chirality2D::write(const MonteCarlo& mc) {
  std::stringstream ss;
  ss << accumulator_->status() << std::endl;
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
