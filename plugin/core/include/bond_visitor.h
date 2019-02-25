
#ifndef FEASST_CORE_BOND_VISITOR_H_
#define FEASST_CORE_BOND_VISITOR_H_

#include <math.h>
#include "core/include/utils_math.h"
#include "core/include/bond.h"
#include "core/include/configuration.h"
#include "core/include/constants.h"

namespace feasst {

class BondTwoBody {
 public:
  virtual double energy(
      const Position& relative,
      const Bond& bond) const = 0;
  virtual ~BondTwoBody() {}
};

/**
  U(r) = 0 when |l-l0| < delta/2, otherwise infinity.
 */
class BondSquareWell : public BondTwoBody {
 public:
  double energy(
      const Position& relative,
      const Bond& bond) const override {
    const double l0 = bond.property("l0");
    const double delta = bond.property("delta");
    if (std::abs(relative.distance() - l0) > 0.5*delta) {
      return NEAR_INFINITY;
    }
    return 0.;
  }
  virtual ~BondSquareWell() {}
};

/**
  angle 0 - 1 - 2
  r01 = r0 - r1, r21 = r2 - r1
 */
class BondThreeBody {
 public:
  virtual double energy(
      const Position& relative01,
      const Position& relative21,
      const Angle& angle) const = 0;
  virtual ~BondThreeBody() {}
};

/**
  U(theta) = 0 when |theta-theta0| < delta/2, otherwise infinity.
  where theta is in degrees
 */
class AngleSquareWell : public BondThreeBody {
 public:
  double energy(
      const Position& relative01,
      const Position& relative21,
      const Angle& angle) const override {
    const double theta0 = angle.property("theta0");
    const double delta = angle.property("delta");
    const double theta = radians_to_degrees(acos(relative01.cosine(relative21)));
    TRACE("theta " << theta);
    if (std::abs(theta - theta0) > 0.5*delta) {
      return NEAR_INFINITY;
    }
    return 0.;
  }
  virtual ~AngleSquareWell() {}
};

class BondVisitor {
 public:
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
        const Angle& angle_type = config.unique_type(part_type).angle(angle.type());
        en += model.energy(relative01, relative21, angle_type);
      }
    }
    set_energy(en);
  }
  void set_energy(const double energy) { energy_ = energy; }
  double energy() const { return energy_; }

 private:
  double energy_;
};

}  // namespace feasst

#endif  // FEASST_CORE_BOND_VISITOR_H_
