#include <cmath>
#include <math.h>  // fabsl
#include "utils/include/serialize.h"
#include "math/include/constants.h"  // PI
#include "math/include/quadratic_equation.h"
#include "math/include/random.h"
#include "math/include/utils_math.h"  // round
#include "monte_carlo/include/trial_select_angle.h"
#include "chain/include/perturb_branch.h"

namespace feasst {

PerturbBranch::PerturbBranch(argtype args) : PerturbBranch(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbBranch::PerturbBranch(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbBranch";
}

class MapPerturbBranch {
 public:
  MapPerturbBranch() {
    auto obj = MakePerturbBranch();
    obj->deserialize_map()["PerturbBranch"] = obj;
  }
};

static MapPerturbBranch mapper_ = MapPerturbBranch();

std::shared_ptr<Perturb> PerturbBranch::create(std::istream& istr) const {
  return std::make_shared<PerturbBranch>(istr);
}

void PerturbBranch::precompute(TrialSelect * select, System * system) {
  ASSERT(system->dimension() == 3, "assumes 3D");
  DEBUG("precomputing a2a1m1");
  a2a1m1_.precompute(select, system);
  DEBUG("precomputing a2a1m2");
  TrialSelectAngle sela2a1m2({{"particle_type", str(select->particle_type())},
    {"mobile_site", feasst::str(select->mobile().site_index(0, 1))},
    {"anchor_site", feasst::str(select->anchor().site_index(0, 0))},
    {"anchor_site2", feasst::str(select->anchor().site_index(0, 1))}});
  sela2a1m2.precompute(system);
  a2a1m2_.precompute(&sela2a1m2, system);
  DEBUG("precomputing m1a1m2");
  TrialSelectAngle selm1a1m2({{"particle_type", str(select->particle_type())},
    {"mobile_site", feasst::str(select->mobile().site_index(0, 0))},
    {"anchor_site", feasst::str(select->anchor().site_index(0, 0))},
    {"anchor_site2", feasst::str(select->mobile().site_index(0, 1))}});
  selm1a1m2.precompute(system);
  m1a1m2_.precompute(&selm1a1m2, system);
}

/*
 * From previous code, translate sites into the following nomenclature:
 * a2 -> 1         1(a2)
 * m1 -> 2         |
 * m2 -> 3         4(a1)  t143(angle)
 * a1 -> 4       /   \L(distance)
 *          2(m1)     3(m2, the site to be placed in this function)
 *
 *               t243(branch_angle)
 *
 * First, check if the branch angles are such that all 4 sites lie on a plane.
 * In such a case, the discriminant of the quadratic below is zero which leads
 * to numerical precision issues.
 * Instead, avoid the quadratic equation in this case.
 *
 * For coordinates of atom 3, x,y,z, letting atom 4 be origin, solve
 * for three eq and 3 unknowns
 * x^2+y^2+z^2 = L^2
 * x*x1+y*y1+z*z1 = cost143   **x1,x2,etc, are unit normal vectors
 * along 14 or 24 bonds
 * x*x2+y*y2+z*z2 = cost243
 * solve for (two values of) x by substitution and quadratic eq. pick
 * one solution randomly
 * if y1 != 0, H = z2 - z1*y2/y1
 * if H != 0, A = (x1y2/y1 - x2)/H, B = (cost243-cost143*y2/y1)/H
 * z(x) = A*x+B, y(x) = C*x+D
 * C = -x1/y1 - Az1/y1, D = cost143/y1 - Bz1/y1
 * [1+C^2+A^2] x^2 + [2CD+2AB] x + [D^2+B^2-L^2] = 0
 *
 * this solution is plagued by numerical stability, if y1 ~ 0 or |H|<1e-8
 *   modified to do alternative solves for more stable H
 * Added a third alternative (thanks Bartosz Mazur).
 */
void PerturbBranch::place_in_branch(const double L,
    const double t143,
    const double t243,
    const double t142,
    System * system,
    TrialSelect * select,
    Random * random) {
  const int dimen = select->configuration(*system).dimension();
  ASSERT(dimen == 3, "assumes 3D, but dimension is " << dimen);

  // check for a planar branch to avoid quadratic equation
  bool planar = false;
  DEBUG("std::abs(t143 + t142 + t243 - 2*PI) " <<
        std::abs(t143 + t142 + t243 - 2*PI) <<
      " std::abs(t143 + t142 - t243) " <<
        std::abs(t143 + t142 - t243));
  if ( (std::abs(t143 + t142 + t243 - 2*PI) < 100*NEAR_ZERO) ||
       (std::abs(t143 + t142 - t243) < 100*NEAR_ZERO) ||
       (std::abs(t243 + t142 - t143) < 100*NEAR_ZERO) ||
       (std::abs(t143 + t243 - t142) < 100*NEAR_ZERO) ) {
    planar = true;
  }

  DEBUG("planar " << planar);
  DEBUG("anchor: " << select->anchor().str());
  DEBUG("anchor site positions size: " << select->anchor().site_positions().size());
  const Position& pos1 = select->anchor_position(0, 1, *system);
  DEBUG("pos1: " << pos1.str());
  const Position& pos4 = select->anchor_position(0, 0, *system);
  DEBUG("pos4: " << pos4.str());
  DEBUG("mobile: " << select->mobile().str());
  const Position& pos2 = select->mobile().site_positions()[0][0];
  DEBUG("pos2: " << pos2.str());
  Position * pos3 = select->get_mobile()->get_site_position(0, 1);
  double x1, y1, z1, x2, y2, z2, r;
  x1 = pos1.coord(0) - pos4.coord(0);
  y1 = pos1.coord(1) - pos4.coord(1);
  z1 = pos1.coord(2) - pos4.coord(2);
  x2 = pos2.coord(0) - pos4.coord(0);
  y2 = pos2.coord(1) - pos4.coord(1);
  z2 = pos2.coord(2) - pos4.coord(2);
  r = std::sqrt(x1*x1+y1*y1+z1*z1);
  x1 /= r; y1 /= r; z1 /= r;
  r = std::sqrt(x2*x2+y2*y2+z2*z2);
  x2 /= r; y2 /= r; z2 /= r;
  const double c143 = std::cos(t143),
               c243 = std::cos(t243);
  double x3, y3, z3;
  DEBUG("x1 " << x1);
  DEBUG("y1 " << y1);
  DEBUG("z1 " << z1);
  DEBUG("H_x " << z2 - x2*z1/x1);
  DEBUG("H_y " << z2 - y2*z1/y1);
  DEBUG("H_z " << x2 - z2*x1/z1);
  const double H_x = z2 - x2*z1/x1;
  const double H_y = z2 - y2*z1/y1;
  const double H_z = x2 - z2*x1/z1;
  bool solved = false;
  double tolerance = 1e-1;
  while (!solved) {
    bool use_h_x = true, use_h_y = true, use_h_z = true;
    if (std::abs(x1) < tolerance) {
      use_h_x = false;
    }
    if (std::abs(y1) < tolerance) {
      use_h_y = false;
    }
    if (std::abs(z1) < tolerance) {
      use_h_z = false;
    }
    if (std::abs(H_x) < tolerance) {
      use_h_x = false;
    }
    if (std::abs(H_y) < tolerance) {
      use_h_y = false;
    }
    if (std::abs(H_z) < tolerance) {
      use_h_z = false;
    }
    if (use_h_y) {
    //if ( std::abs(H_y) >= std::abs(H_x) && std::abs(H_y) >= std::abs(H_z) ) {
      DEBUG("use H_y");
    //if ( std::abs(y1) > std::abs(x1) ) {
  //  if ( (abs(x1) < DTOL) || (abs(Cyz) > abs(Cxz)) ) {
      // cout << "test1 " << abs(x1) << " t2 " << abs(Cyz) << " > "
      //      << abs(Cxz) << endl;
      solve_branch_(x1, y1, z1, x2, y2, z2, &x3, &y3, &z3, c143, c243, planar, random);
      solved = true;
    } else if (use_h_x) {
    //} else if (std::abs(H_x) >= std::abs(H_z)) {
      DEBUG("use H_x");
    // } else if ( (abs(y1) < DTOL) || (abs(Cyz) < abs(Cxz)) ) {
      // cout << "test2 " << abs(y1) << " t2 " << abs(Cyz) << " > "
      //      << abs(Cxz) << endl;
      solve_branch_(y1, x1, z1, y2, x2, z2, &y3, &x3, &z3, c143, c243, planar, random);
      solved = true;
    } else if (use_h_z) {
    //} else {
      DEBUG("use H_z");
      solve_branch_(y1, z1, x1, y2, z2, x2, &y3, &z3, &x3, c143, c243, planar, random);
      solved = true;
    }
    tolerance /= 10.;
    if (tolerance < 1e-15) {
      FATAL("PerturbBranch error. " << x1 << " " << y1 << " " << z1 << " " <<
        x2 << " " << y2 << " " << z2 << " " << c143 << " " << c243 << " " <<
        planar);
    }
  }
  //DEBUG("L: " << L << " dist " << std::sqrt(x3*x3+y3*y3+z3*z3));
//  const double dist = std::sqrt(x3*x3+y3*y3+z3*z3);
//  x3 /= dist;
//  y3 /= dist;
//  z3 /= dist;
  pos3->set_coord(0, L*x3 + pos4.coord(0));
  pos3->set_coord(1, L*y3 + pos4.coord(1));
  pos3->set_coord(2, L*z3 + pos4.coord(2));
  select->get_configuration(system)->update_positions(select->mobile());
}

void PerturbBranch::solve_branch_(
    const double x1, const double y1, const double z1,
    const double x2, const double y2, const double z2, double *x3, double *y3,
    double *z3, const double c143, const double c243, const bool planar,
    Random * random) const {
  ASSERT(y1 != 0, "y1==0");
  const long double H = z2 - y2*z1/y1;
  ASSERT(H != 0, "H==0");
//  DEBUG("H " << H);
  const long double A = (x1*y2/y1 - x2)/H,
               B = (c243 - c143*y2/y1)/H,
               C = -x1/y1 - A*z1/y1,
               D = c143/y1 - B*z1/y1,
               a = (1+A*A+C*C),
               b = 2*(A*B+C*D),
               c = (B*B+D*D-1);
  long double ans1 = NEAR_INFINITY;
//  DEBUG("a " << a);
//  DEBUG("b " << b);
//  DEBUG("c " << c);
  if (planar) {
    ans1 = -b/2/a;
  } else {
    long double ans2 = NEAR_INFINITY;
    long double discrim = -1;
    quadratic_equation(a, b, c, &discrim, &ans1, &ans2);
  //  DEBUG("discrim " << discrim);
  //  DEBUG("ans1 " << ans1);
  //  DEBUG("ans2 " << ans2);
    if (discrim < 0) {
  //    DEBUG("discrim " << discrim);
      if ( (std::sqrt(fabsl(discrim))/2/fabsl(a) < 1000*std::sqrt(NEAR_ZERO)) ||
           fabsl(discrim) < 10000*NEAR_ZERO) {
        // within double preicison, the discriminant is zero
        ans1 = ans2 = -b/2/a;
      } else {
        std::streamsize ss = std::cout.precision();
        std::cout << std::setprecision(std::numeric_limits<long double>::digits10+2)
             << "c143 " << c143 << " c243 " << c243 << std::endl;
        std::cout << "x1 " << x1 << " " << y1 << " " << z1 << std::endl;
        std::cout << "x2 " << x2 << " " << y2 << " " << z2 << std::endl;
        std::cout << "A " << A << " H " << H << " B " << B << " C " << C
             << " D " << D << std::endl;
        std::cout << "ans1 " << ans1 << " ans2 " << ans2 << std::endl;
        std::cout << "discrim " << discrim << std::endl;
        std::cout << "a " << a << " b " << b << " c " << c << std::endl;
        std::cout << std::setprecision(ss);
        std::cout << "tol " << std::sqrt(NEAR_ZERO) << " rel "
             << std::sqrt(fabsl(discrim))/2/fabsl(a) << std::endl;
        std::cout << "tol " << 10*NEAR_ZERO << " rel " << fabsl(discrim)
             << std::endl;
        FATAL("imaginary branch");
      }
    }
    if (random->coin_flip()) {
      ans1 = ans2;
    }
  }

  // (x,y,z) is unit vector pointing in direction of l3 bond
  *x3 = ans1;
  *y3 = C*ans1 + D;
  *z3 = A*ans1 + B;
}

void PerturbBranch::move(const bool is_position_held,
                         System * system,
                         TrialSelect * select,
                         Random * random) {
  const Angle& a2a1m1 = select->configuration(*system).unique_type(
    select->particle_type()).angle(a2a1m1_.angle_type());
  const Angle& a2a1m2 = select->configuration(*system).unique_type(
    select->particle_type()).angle(a2a1m2_.angle_type());
  const Angle& m1a1m2 = select->configuration(*system).unique_type(
    select->particle_type()).angle(m1a1m2_.angle_type());
  double radians_a2a1m1, radians_a2a1m2, radians_m1a1m2;
  ASSERT(angle_.deserialize_map().count(a2a1m1.model()) == 1,
    a2a1m1.model() << " not found");
  const BondThreeBody * model = angle_.deserialize_map()[a2a1m1.model()].get();
  double bond_energy = 0.;
  if (is_position_held) {
    const Position& a1 = select->anchor_position(0, 0, *system);
    const Position& a2 = select->anchor_position(0, 1, *system);
    const Position& m1 = select->mobile().site_positions()[0][0];
    const Position& m2 = select->mobile().site_positions()[0][1];
    bond_energy += model->energy(a2, a1, m1, a2a1m1);
    bond_energy += model->energy(a2, a1, m2, a2a1m2);
    bond_energy += model->energy(m1, a1, m2, m1a1m2);
  } else {
    model->random_branch(
      a2a1m1, a2a1m2, m1a1m2,
      system->thermo_params().beta(),
      &radians_a2a1m1, &radians_a2a1m2, &radians_m1a1m2,
      random);
    bond_energy += model->energy(radians_a2a1m1, a2a1m1);
    bond_energy += model->energy(radians_a2a1m2, a2a1m2);
    bond_energy += model->energy(radians_m1a1m2, m1a1m2);
    const double la1m1 = a2a1m1_.random_distance(*system, select, random, &bond_energy);
//    DEBUG("la1m1 " << la1m1);
    const double la1m2 = a2a1m2_.random_distance(*system, select, random, &bond_energy);
//    DEBUG("la1m2 " << la1m2);
    a2a1m1_.place_in_circle(la1m1, radians_a2a1m1, system, select, random);
    place_in_branch(la1m2, radians_a2a1m2, radians_m1a1m2, radians_a2a1m1, system, select, random);
//    DEBUG("here");
//    for (const std::vector<Position>& poss : select->anchor().site_positions()) {
//      for (const Position& pos : poss) {
//        DEBUG(pos.str());
//      }
//    }
//    for (const std::vector<Position>& poss : select->mobile().site_positions()) {
//      for (const Position& pos : poss) {
//        DEBUG(pos.str());
//      }
//    }
  }
  select->add_exclude_energy(bond_energy);
}

PerturbBranch::PerturbBranch(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbBranch", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(3905 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&a2a1m1_, istr);
  feasst_deserialize_fstobj(&a2a1m2_, istr);
  feasst_deserialize_fstobj(&m1a1m2_, istr);
}

void PerturbBranch::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(3905, ostr);
  feasst_serialize_fstobj(a2a1m1_, ostr);
  feasst_serialize_fstobj(a2a1m2_, ostr);
  feasst_serialize_fstobj(m1a1m2_, ostr);
}

}  // namespace feasst
