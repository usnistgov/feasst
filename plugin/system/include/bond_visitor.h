
#ifndef FEASST_CONFIGURATION_BOND_VISITOR_H_
#define FEASST_CONFIGURATION_BOND_VISITOR_H_

#include <math.h>
#include <memory>
#include <string>
#include <map>
#include "math/include/utils_math.h"
#include "configuration/include/bond.h"
#include "configuration/include/configuration.h"
#include "system/include/bond_two_body.h"
#include "system/include/bond_three_body.h"

namespace feasst {

class BondVisitor {
 public:
  explicit BondVisitor(const argtype& args = argtype()) {}
  void compute(
      const BondTwoBody& model,
      const Configuration& config,
      const int group_index = 0) {
    const Select& selection = config.group_selects()[group_index];
    compute(model, selection, config);
  }
  void compute(
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
  void compute(
      const BondThreeBody& model,
      const Configuration& config,
      const int group_index = 0) {
    const Select& selection = config.group_selects()[group_index];
    compute(model, selection, config);
  }
  void compute(
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
  void set_energy(const double energy) { energy_ = energy; }
  double energy() const { return energy_; }

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<BondVisitor> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<BondVisitor> >& deserialize_map();
  std::shared_ptr<BondVisitor> deserialize(std::istream& istr);
  explicit BondVisitor(std::istream& istr);
  virtual ~BondVisitor() {}

 protected:
  std::string class_name_ = "BondVisitor";

  void serialize_bond_visitor_(std::ostream& ostr) const;

 private:
  double energy_;
};

inline std::shared_ptr<BondVisitor> MakeBondVisitor(
    const argtype &args = argtype()) {
  return std::make_shared<BondVisitor>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_VISITOR_H_
