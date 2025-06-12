#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/position.h"
#include "math/include/utils_math.h"
#include "configuration/include/bond.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "chain/include/analyze_bonds.h"

namespace feasst {

FEASST_MAPPER(AnalyzeBonds,);

AnalyzeBonds::AnalyzeBonds(argtype * args) : Analyze(args) {
  Histogram bhist, ahist, dhist;
  bhist.set_width_center(dble("bond_bin_width", args, 1),
                         dble("bond_bin_center", args, 0.));
  bond_hist_.push_back(bhist);
  ahist.set_width_center(dble("angle_bin_width", args, 0.01),
                         dble("angle_bin_center", args, 0.));
  angle_hist_.push_back(ahist);
  dhist.set_width_center(dble("dihedral_bin_width", args, 0.01),
                         dble("dihedral_bin_center", args, 0.));
  dihedral_hist_.push_back(dhist);
}
AnalyzeBonds::AnalyzeBonds(argtype args) : AnalyzeBonds(&args) {
  feasst_check_all_used(args);
}

void AnalyzeBonds::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2390, ostr);
  feasst_serialize_fstobj(bond_, ostr);
  feasst_serialize_fstobj(angle_, ostr);
  feasst_serialize_fstobj(dihedral_, ostr);
  feasst_serialize_fstobj(bond_hist_, ostr);
  feasst_serialize_fstobj(angle_hist_, ostr);
  feasst_serialize_fstobj(dihedral_hist_, ostr);
}

AnalyzeBonds::AnalyzeBonds(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2390, "version mismatch: " << version);
  feasst_deserialize_fstobj(&bond_, istr);
  feasst_deserialize_fstobj(&angle_, istr);
  feasst_deserialize_fstobj(&dihedral_, istr);
  feasst_deserialize_fstobj(&bond_hist_, istr);
  feasst_deserialize_fstobj(&angle_hist_, istr);
  feasst_deserialize_fstobj(&dihedral_hist_, istr);
}

void AnalyzeBonds::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  DEBUG("initializing AnalyzeBonds");
  const Configuration& config = mc->system().configuration();
  for (int btype = 0; btype < config.num_bond_types(); ++btype) {
    if (static_cast<int>(bond_.size()) != config.num_bond_types()) {
      bond_.push_back(Accumulator());
      if (btype != 0) bond_hist_.push_back(bond_hist_[0]);
    }
  }
  for (int atype = 0; atype < config.num_angle_types(); ++atype) {
    if (static_cast<int>(angle_.size()) != config.num_angle_types()) {
      angle_.push_back(Accumulator());
      if (atype != 0) angle_hist_.push_back(angle_hist_[0]);
    }
  }
  for (int dtype = 0; dtype < config.num_dihedral_types(); ++dtype) {
    if (static_cast<int>(dihedral_.size()) != config.num_dihedral_types()) {
      dihedral_.push_back(Accumulator());
      if (dtype != 0) dihedral_hist_.push_back(dihedral_hist_[0]);
    }
  }
}

void AnalyzeBonds::update(const MonteCarlo& mc) {
  const Configuration& config = mc.system().configuration();
  const Select& all = config.selection_of_all();
  for (int index = 0; index < all.num_particles(); ++index) {
    const int part_index = all.particle_index(index);
    const Particle& part = config.select_particle(part_index);
    const int part_type = part.type();
    for (const Bond& bond : config.particle_type(part_type).bonds()) {
      const Position& ri = part.site(bond.site(0)).position();
      const Position& rj = part.site(bond.site(1)).position();
      const double distance = ri.distance(rj);
      bond_[bond.type()].accumulate(distance);
      bond_hist_[bond.type()].add(distance);
    }
    for (const Angle& angle : config.particle_type(part_type).angles()) {
      const Position& ri = part.site(angle.site(0)).position();
      const Position& rj = part.site(angle.site(1)).position();
      const Position& rk = part.site(angle.site(2)).position();
      const double theta = rj.vertex_angle_radians(ri, rk);
      angle_[angle.type()].accumulate(theta);
      angle_hist_[angle.type()].add(theta);
    }
    for (const Dihedral& dihedral : config.particle_type(part_type).dihedrals()) {
      const Position& ri = part.site(dihedral.site(0)).position();
      const Position& rj = part.site(dihedral.site(1)).position();
      const Position& rk = part.site(dihedral.site(2)).position();
      const Position& rl = part.site(dihedral.site(3)).position();
      const double phi = ri.torsion_angle_radians(rj, rk, rl);
      dihedral_[dihedral.type()].accumulate(phi);
      dihedral_hist_[dihedral.type()].add(phi);
    }
  }
}

std::string AnalyzeBonds::write(const MonteCarlo& mc) {
  std::stringstream ss;
  for (int type = 0; type < static_cast<int>(bond_.size()); ++type) {
    DEBUG("type " << type);
    DEBUG("bond size " << bond_.size());
    ss << "# bond: " << type << ": "  << bond_[type].str() << std::endl;
    DEBUG("bond hist size " << bond_hist_.size());
    ss << "# bond hist: " << type << ": " << bond_hist_[type].str() << std::endl;
  }
  for (int type = 0; type < static_cast<int>(angle_.size()); ++type) {
    ss << "# angle: " << type << ": "  << angle_[type].str() << std::endl;
    ss << "# angle hist: " << type << ": " << angle_hist_[type].str() << std::endl;
  }
  for (int type = 0; type < static_cast<int>(dihedral_.size()); ++type) {
    ss << "# dihedral: " << type << ": "  << dihedral_[type].str() << std::endl;
    ss << "# dihedral hist: " << type << ": " << dihedral_hist_[type].str() << std::endl;
  }
  return ss.str();
}

}  // namespace feasst
