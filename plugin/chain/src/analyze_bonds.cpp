#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "chain/include/analyze_bonds.h"

namespace feasst {

class MapAnalyzeBonds {
 public:
  MapAnalyzeBonds() {
    AnalyzeBonds().deserialize_map()["AnalyzeBonds"] =
      MakeAnalyzeBonds();
  }
};

static MapAnalyzeBonds mapper_ = MapAnalyzeBonds();

AnalyzeBonds::AnalyzeBonds(argtype * args) : AnalyzeUpdateOnly(args) {
  Histogram bhist, ahist, dhist;
  bhist.set_width_center(dble("bond_bin_width", args, 1),
                         dble("bond_bin_center", args, 0.));
  bond_hist_.push_back(bhist);
  ahist.set_width_center(dble("angle_bin_width", args, 1),
                         dble("angle_bin_center", args, 0.));
  angle_hist_.push_back(ahist);
  dhist.set_width_center(dble("dihedral_bin_width", args, 1),
                         dble("dihedral_bin_center", args, 0.));
  dihedral_hist_.push_back(dhist);
}
AnalyzeBonds::AnalyzeBonds(argtype args) : AnalyzeBonds(&args) {
  FEASST_CHECK_ALL_USED(args);
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

AnalyzeBonds::AnalyzeBonds(std::istream& istr)
  : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2390, "version mismatch: " << version);
  feasst_deserialize_fstobj(&bond_, istr);
  feasst_deserialize_fstobj(&angle_, istr);
  feasst_deserialize_fstobj(&dihedral_, istr);
  feasst_deserialize_fstobj(&bond_hist_, istr);
  feasst_deserialize_fstobj(&angle_hist_, istr);
  feasst_deserialize_fstobj(&dihedral_hist_, istr);
}

void AnalyzeBonds::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  const Configuration& config = system->configuration();
  for (int btype = 0; btype < config.num_bond_types(); ++btype) {
    bond_.push_back(Accumulator());
    if (btype != 0) bond_hist_.push_back(bond_hist_[0]);
  }
  for (int atype = 0; atype < config.num_angle_types(); ++atype) {
    angle_.push_back(Accumulator());
    if (atype != 0) angle_hist_.push_back(angle_hist_[0]);
  }
  for (int dtype = 0; dtype < config.num_dihedral_types(); ++dtype) {
    dihedral_.push_back(Accumulator());
    if (dtype != 0) dihedral_hist_.push_back(dihedral_hist_[0]);
  }
}

void AnalyzeBonds::update(const Criteria& criteria,
  const System& system,
  const TrialFactory& trial_factory) {
  const Configuration& config = system.configuration();
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

}  // namespace feasst
