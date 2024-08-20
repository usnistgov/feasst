#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h" // round
#include "math/include/constants.h"
#include "math/include/random.h" // normal
#include "math/include/position.h"
#include "configuration/include/bond.h"
#include "models/include/angle_harmonic.h"

namespace feasst {

class MapAngleHarmonic {
 public:
  MapAngleHarmonic() {
    auto obj = MakeAngleHarmonic();
    obj->deserialize_map()["AngleHarmonic"] = obj;
  }
};

static MapAngleHarmonic mapper_ = MapAngleHarmonic();

std::shared_ptr<BondThreeBody> AngleHarmonic::create(std::istream& istr) const {
  return std::make_shared<AngleHarmonic>(istr);
}

AngleHarmonic::AngleHarmonic(std::istream& istr) : BondThreeBody(istr) {
  // ASSERT(class_name_ == "AngleHarmonic", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(4655 == version, "mismatch version: " << version);
}

void AngleHarmonic::serialize_angle_harmonic_(std::ostream& ostr) const {
  serialize_bond_three_body_(ostr);
  feasst_serialize_version(4655, ostr);
}

void AngleHarmonic::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_angle_harmonic_(ostr);
}

double AngleHarmonic::energy(const double radians, const Bond& angle) const {
  DEBUG("radians " << radians);
  const double equil_radians =
    degrees_to_radians(angle.property("equilibrium_degrees"));
  const double k = angle.property("k_energy_per_radian_sq");
  double delta_rad = radians - equil_radians;
  DEBUG("delta_rad " << delta_rad);
  return k*delta_rad*delta_rad;
}

//double AngleHarmonic::random_angle_radians(const Angle& angle, const double beta,
//    const int dimension, Random * random) const {
//  if (dimension == 2) {
//    return BondThreeBody::random_angle_radians(angle, beta, dimension, random);
//  }
//  ASSERT(dimension == 3, "unrecognized dimension: " << dimension);
//  double min_rad = 0.;
//  if (angle.has_property("minimum_degrees")) {
//    min_rad = degrees_to_radians(angle.property("minimum_degrees"));
//  }
//  const double equil_radians =
//    degrees_to_radians(angle.property("equilibrium_degrees"));
//  const double k = angle.property("k_energy_per_radian_sq");
//  int attempt = 0;
//  while (attempt < 1e6) {
//    const double radians = random->normal(equil_radians,
//      std::sqrt(beta*k*0.5));
//    if (radians >= min_rad && radians <= PI) {
//      const double en = energy(radians, angle);
//      const double jacobian = std::sin(radians);
//      if (random->uniform() < jacobian*std::exp(-beta*en)) {
//        return radians;
//      }
//    }
//    ++attempt;
//  }
//  FATAL("max attempts reached");
//}

// Implementation of two-branch num_jacobian_gaussian
// attempt a number of branches
// select based on jacobian-gaussian weight
// Note that torsion angle reflection, theta -> 2pi-theta
//   has invariant |sin(theta)| because sin(2pi - theta) = -sin(theta)
//   simply compute one of the torsion angles for jacobian
// add rosenbluth to acceptance via ln_met

void AngleHarmonic::random_branch(
    const Angle& a2a1m1,
    const Angle& a2a1m2,
    const Angle& m1a1m2,
    const double beta,
    const bool is_position_held,
    double * radians_a2a1m1,
    double * radians_a2a1m2,
    double * radians_m1a1m2,
    Random * random,
    double * ln_met,
    const Position * const a1,
    const Position * const a2,
    const Position * const m1,
    const Position * const m2) const {
  DEBUG("is_position_held " << is_position_held);
  // See PerturbBranch for a1, a2, m1, m2 definition
  const double theta1_eq = degrees_to_radians(a2a1m1.property("equilibrium_degrees"));
  const double theta2_eq = degrees_to_radians(a2a1m2.property("equilibrium_degrees"));
  const double theta12_eq = degrees_to_radians(m1a1m2.property("equilibrium_degrees"));
  const double k1 = a2a1m1.property("k_energy_per_radian_sq");
  const double k2 = a2a1m2.property("k_energy_per_radian_sq");
  const double k12 = m1a1m2.property("k_energy_per_radian_sq");
  int num_jacobian_gaussian = 0;
  if (m1a1m2.has_property("num_jacobian_gaussian")) {
    num_jacobian_gaussian = round(m1a1m2.property("num_jacobian_gaussian"));
  }
  DEBUG("num_jacobian_gaussian " << num_jacobian_gaussian);
  if (num_jacobian_gaussian <= 0) {
    BondThreeBody::random_branch(a2a1m1, a2a1m2, m1a1m2, beta, is_position_held,
                                 radians_a2a1m1, radians_a2a1m2, radians_m1a1m2,
                                 random, ln_met, a1, a2, m1, m2);
    //WARN("not using Jacobian-Gaussian");
    return;
  }
  struct jg {
    double t1, t2, t12;
  };
  std::vector<jg> steps(num_jacobian_gaussian);
  std::vector<double> weight(num_jacobian_gaussian);
  //for (jg& step : steps) {
  for (int istep = 0; istep < num_jacobian_gaussian; ++istep) {
    jg& step = steps[istep];
    if (istep == 0 && is_position_held) {
      ASSERT(a1, "error");
      step.t1 = a1->vertex_angle_radians(*a2, *m1);
      step.t2 = a1->vertex_angle_radians(*a2, *m2);
      step.t12 = a1->vertex_angle_radians(*m1, *m2);
    } else {
      int attempt = 0;
      bool accept = false;
      while (!accept) {
        if (std::abs(k1) < NEAR_ZERO) {
          step.t1 = PI*random->uniform();
        } else {
          step.t1 = random->normal(theta1_eq, 1./std::sqrt(beta*2.*k1));
        }
        DEBUG("trial t1 " << step.t1);
        if (step.t1 >= 0. && step.t1 <= PI) {
          if (std::abs(k2) < NEAR_ZERO) {
            step.t2 = PI*random->uniform();
          } else {
            step.t2 = random->normal(theta2_eq, 1./std::sqrt(beta*2.*k2));
          }
          DEBUG("trial t2 " << step.t2);
          if (step.t2 >= 0. && step.t2 <= PI) {
            if (std::abs(k12) < NEAR_ZERO) {
              step.t12 = PI*random->uniform();
            } else {
              step.t12 = random->normal(theta12_eq, 1./std::sqrt(beta*2.*k12));
            }
            DEBUG("trial t12 " << step.t12);
            if (step.t12 >= std::abs(step.t1 - step.t2) &&
                step.t12 <= std::min(step.t1 + step.t2, 2.*PI-(step.t1 + step.t2))) {
              accept = true;
            }
          }
        }
        ++attempt;
        if (attempt == 1e6) FATAL("max attempts reached");
      }
    }
    DEBUG("t1 " << step.t1);
    DEBUG("t2 " << step.t2);
    DEBUG("t12 " << step.t12);
    double w12 = std::cos(step.t12) - std::cos(step.t1)*std::cos(step.t2);
    w12 = std::acos(w12/std::sin(step.t1)/std::sin(step.t2));
    DEBUG("w12 " << w12);
    // w12 may be randomly reflected, w12->2pi-w12
    // however, the |sin(w12)| term in jacobian is unaffected
    // because sin(2pi-t)=-sin(t)
    weight[istep] = std::abs(std::sin(step.t12)/std::sin(w12));
    if (istep > 0) weight[istep] += weight[istep - 1];
  }
  const double sum = weight.back();
  // met prob seems to be missing a factor of 0.82 based on GCE tests
  *ln_met = std::log(sum/static_cast<double>(num_jacobian_gaussian));
  if (is_position_held) {
    *ln_met *= -1;
  };
  for (double& w : weight) w /= weight.back();
//  DEBUG(feasst_str(weight));
  const int index = random->index_from_cumulative_probability(weight);
  *radians_a2a1m1 = steps[index].t1;
  *radians_a2a1m2 = steps[index].t2;
  *radians_m1a1m2 = steps[index].t12;
}

}  // namespace feasst
