
#ifndef FEASST_CONFIGURATION_BOND_THREE_BODY_H_
#define FEASST_CONFIGURATION_BOND_THREE_BODY_H_

#include <map>
#include <string>
#include <memory>
#include "math/include/position.h"
#include "configuration/include/bond.h"

namespace feasst {

class Random;

/**
  A three body bond is defined by three sites, 0 - 1 - 2.
  The relative vector r01 = r0 - r1 points from r1 to r0.
  The relative vector r21 = r2 - r1 points from r1 to r2.

  Angle Property minimum_degrees sets the energy to NEAR_INFINITY at smaller
  angles.
 */
class BondThreeBody {
 public:
  BondThreeBody() {}
  virtual double energy(const Position& ri, const Position& rj,
    const Position& rk, const Bond& angle) const;
  virtual double energy(const double radians, const Bond& angle) const = 0;

  /**
    Return a randomly selected bond angle.

    \f$P(\theta) \propto \sin\theta\exp[-\beta U(\theta)]\f$
   */
  virtual double random_angle_radians(const Angle& angle, const double beta,
    const int dimension, Random * random) const;

  /**
    Return three random angles for forming a branch.

    anchor2 -> 1         1(a2)
    mobile1 -> 2         |
    mobile2 -> 3         4(a1)  t143(angle)
    anchor1 -> 4       /   \L34(distance)
                  2(m1)     3(m2, the site to be placed in this function)
                       t243(branch_angle)
   */
  virtual void random_branch(
    const Angle& a2a1m1,
    const Angle& a2a1m2,
    const Angle& m1a1m2,
    const double beta,
    double * radians_a2a1m1,
    double * radians_a2a1m2,
    double * radians_m1a1m2,
    Random * random) const;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<BondThreeBody> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<BondThreeBody> >& deserialize_map();
  std::shared_ptr<BondThreeBody> deserialize(std::istream& istr);
  virtual ~BondThreeBody() {}

 protected:
  std::string class_name_ = "BondThreeBody";

  void serialize_bond_three_body_(std::ostream& ostr) const;
  explicit BondThreeBody(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_THREE_BODY_H_
