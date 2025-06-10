
#ifndef FEASST_SUPERTOROID_SPHERE_H_
#define FEASST_SUPERTOROID_SPHERE_H_

#include "shape/include/shape.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Random;

/**
  As currently implemented, the is_inside function does not take into account
  the diameter of the particle.
  So different size particles will experience difference sized cavities if
  using ModelHardShape.
  Put another way, only the center position is considered.

  A supertoroid is given by the following surface vector

  \f$\left\{
    \begin{array}{lr}
      a_1[a_4+\cos(\eta)^\epsilon_1]\cos(\omega)^\epsilon_2 & -\pi \le \eta \le \pi \\
      a_2[a_4+\cos(\eta)^\epsilon_1]\sin(\omega)^\epsilon_2 & -\pi \le \omega \le \pi \\
      a_3\sin(\eta)^\epsilon_1 &
    \end{array}
  \right\}\f$

  with an implicit function, \f$F\f$, where \f$F=1\f$ at the surface,
  \f$F<1\f$ inside, and \f$F>1\f$ outside,

  \f$\left(\left(\left(\frac{x}{a_1}\right)^{\frac{2}{\epsilon_2}} + \left(\frac{y}{a_2}\right)^{\frac{2}{\epsilon_2}}\right)^{\frac{\epsilon_2}{2}}-a_4\right)^{\frac{2}{\epsilon_1}}+\left(\frac{z}{a_3}\right)^{\frac{2}{\epsilon_1}}\f$.

  The size of the hole is related to \f$a_4\f$, where the radius of the toroid, \f$R=a_4\sqrt{a_1^2+a_2^2}\f$.

  To model a superquadric, set \f$a_4=0\f$.
 */
class Supertoroid : public Shape {
 public:
  //@{
  /** @name Arguments
    - center: comma-separated values for the positions in each dimension.
      (default: a three dimensional origin is assumed, e.g., center=0,0,0).
    - a1: (default: 1).
    - a2: (default: 1).
    - a3: (default: 1).
    - a4: (default: 0).
    - epsilon1: (default: 1).
    - epsilon2: (default: 1).
   */
  explicit Supertoroid(argtype args = argtype());
  explicit Supertoroid(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  const Position& center() const { return center_; }
  double nearest_distance(const Position& point) const override;
  bool is_inside(const Position& point) const override;
  bool is_inside(const Position& point, const double diameter) const override;
  double surface_area() const override;
  double volume() const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<Supertoroid>(istr); }
  std::shared_ptr<Shape> create(argtype * args) const override {
    return std::make_shared<Supertoroid>(args); }
  explicit Supertoroid(std::istream& istr);
  virtual ~Supertoroid() {}

  //@}
 private:
  Position center_;
  double a1_, a2_, a3_, a4_, epsilon1_, epsilon2_;
};

inline std::shared_ptr<Supertoroid> MakeSupertoroid(argtype args = argtype()) {
  return std::make_shared<Supertoroid>(args);
}

}  // namespace feasst

#endif  // FEASST_SUPERTOROID_SPHERE_H_
