#include "utils/include/serialize.h"
#include "configuration/include/configuration.h"
#include "system/include/bond_visitor.h"

namespace feasst {

class MapBondVisitor {
 public:
  MapBondVisitor() {
    auto obj = MakeBondVisitor();
    obj->deserialize_map()["BondVisitor"] = obj;
  }
};

static MapBondVisitor mapper_ = MapBondVisitor();

std::map<std::string, std::shared_ptr<BondVisitor> >& BondVisitor::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondVisitor> >* ans =
     new std::map<std::string, std::shared_ptr<BondVisitor> >();
  return *ans;
}

void BondVisitor::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bond_visitor_(ostr);
}

std::shared_ptr<BondVisitor> BondVisitor::create(std::istream& istr) const {
  return std::make_shared<BondVisitor>(istr);
}

std::shared_ptr<BondVisitor> BondVisitor::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondVisitor::serialize_bond_visitor_(std::ostream& ostr) const {
  feasst_serialize_version(303, ostr);
  feasst_serialize(energy_, ostr);
}

BondVisitor::BondVisitor(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(303 == version, "mismatch version: " << version);
  feasst_deserialize(&energy_, istr);
}

void BondVisitor::compute(
    const BondTwoBody& model,
    const Select& selection,
    const Configuration& config) {
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const int part_type = part.type();
    for (const Bond& bond : config.particle_type(part_type).bonds()) {
      const Position& position0 = part.site(bond.site(0)).position();
      const Position& position1 = part.site(bond.site(1)).position();
      Position relative = position0;
      relative.subtract(position1);
      TRACE("bond ij " << part_index << " " << bond.site(0) << " "
        << bond.site(1) << " sq " << relative.squared_distance());
      const Bond& bond_type = config.unique_type(part_type).bond(bond.type());
      en += model.energy(relative, bond_type);
    }
  }
  set_energy(en);
}

void BondVisitor::compute(
    const BondThreeBody& model,
    const Select& selection,
    const Configuration& config) {
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const int part_type = part.type();
    for (const Angle& angle : config.particle_type(part_type).angles()) {
      const Position& position0 = part.site(angle.site(0)).position();
      const Position& position1 = part.site(angle.site(1)).position();
      const Position& position2 = part.site(angle.site(2)).position();
      Position relative01 = position0;
      Position relative21 = position2;
      relative01.subtract(position1);
      relative21.subtract(position1);
      const Angle& angle_type =
        config.unique_type(part_type).angle(angle.type());
      en += model.energy(relative01, relative21, angle_type);
    }
  }
  set_energy(en);
}

void BondVisitor::compute(
    const BondTwoBody& model,
    const Configuration& config,
    const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute(model, selection, config);
}

void BondVisitor::compute(
    const BondThreeBody& model,
    const Configuration& config,
    const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute(model, selection, config);
}
}  // namespace feasst
