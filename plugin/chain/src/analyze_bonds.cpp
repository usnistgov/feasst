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

AnalyzeBonds::AnalyzeBonds(const argtype &args)
  : AnalyzeUpdateOnly(args) {
  Arguments args_(args);
  args_.dont_check();
  Histogram bhist, ahist;
  bhist.set_width_center(args_.key("bond_bin_width").dflt("1").dble(), 0.);
  bond_hist_.push_back(bhist);
  ahist.set_width_center(args_.key("angle_bin_width").dflt("1").dble(), 0.);
  angle_hist_.push_back(ahist);
}

void AnalyzeBonds::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2390, ostr);
  feasst_serialize_fstobj(bond_, ostr);
  feasst_serialize_fstobj(angle_, ostr);
  feasst_serialize_fstobj(bond_hist_, ostr);
  feasst_serialize_fstobj(angle_hist_, ostr);
}

AnalyzeBonds::AnalyzeBonds(std::istream& istr)
  : AnalyzeUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2390, "version mismatch: " << version);
  feasst_deserialize_fstobj(&bond_, istr);
  feasst_deserialize_fstobj(&angle_, istr);
  feasst_deserialize_fstobj(&bond_hist_, istr);
  feasst_deserialize_fstobj(&angle_hist_, istr);
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
      const Position& position0 = part.site(bond.site(0)).position();
      const Position& position1 = part.site(bond.site(1)).position();
      Position relative = position0;
      relative.subtract(position1);
      const double distance = relative.distance();
      bond_[bond.type()].accumulate(distance);
      bond_hist_[bond.type()].add(distance);
    }
    for (const Angle& angle : config.particle_type(part_type).angles()) {
      const Position& position0 = part.site(angle.site(0)).position();
      const Position& position1 = part.site(angle.site(1)).position();
      const Position& position2 = part.site(angle.site(2)).position();
      Position relative01 = position0;
      Position relative21 = position2;
      relative01.subtract(position1);
      relative21.subtract(position1);
      const double theta = radians_to_degrees(
        std::acos(relative01.cosine(relative21)));
      angle_[angle.type()].accumulate(theta);
      angle_hist_[angle.type()].add(theta);
    }
  }
}

}  // namespace feasst
