#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/formula.h"
#include "math/include/constants.h"
// #include "math/include/solver_newton_raphson.h"
// #include "math/include/solver_bisection.h"
#include "math/include/solver_brent_dekker.h"
#include "shape/include/half_space_sine.h"

namespace feasst {

class MapHalfSpaceSine {
 public:
  MapHalfSpaceSine() {
    auto obj = MakeHalfSpaceSine(
      MakeFormulaSineWave(),
      { {"dimension", "1"},
        {"intersection", "1"},
        {"direction", "1"},
        {"wave_dimension", "0"}});
    obj->deserialize_map()["HalfSpaceSine"] = obj;
  }
};

static MapHalfSpaceSine mapper_ = MapHalfSpaceSine();

// https://math.stackexchange.com/questions/514820/distance-between-point-and-sine-wave#:~:text=The%20points%20on%20the%20sine,its%20derivative%20will%20be%20zero.
// HWH.. root finder for nearest distance?
// HWH.. how about not use nearest_distance.. just use is_inside?
class SineDistDeriv : public Formula {
 public:
  SineDistDeriv(const HalfSpaceSine& wall) {
    wall_ = const_cast<const HalfSpaceSine *>(&wall);
    wave_ = const_cast<const FormulaSineWave *>(&wall.sine_wave());
  }
  void set_point(const Position& point) {
    point_ = const_cast<const Position *>(&point);
  }

  double evaluate(const double x) const override {
    const double e = point_->coord(wall_->wave_dimension());
    const double f = point_->coord(wall_->dimension());
    return x - e + (wave_->evaluate(x) - f)*(wave_->derivative(x));
  }
  double derivative(const double x) const override {
    const double f = point_->coord(wall_->dimension());
    return 1. +wave_->second_derivative(x)*(wave_->evaluate(x) - f)
         + std::pow(wave_->derivative(x), 2);
  }
 private:
  const FormulaSineWave * wave_;
  const HalfSpaceSine * wall_;
  const Position * point_;
};

void HalfSpaceSine::init_sine_dist_deriv_() {
  sine_dist_deriv_ = std::make_shared<SineDistDeriv>(*this);
}

HalfSpaceSine::HalfSpaceSine(std::shared_ptr<FormulaSineWave> sine_wave,
    argtype args) : HalfSpace(&args) {
  class_name_ = "HalfSpaceSine";
  wave_dimension_ = integer("wave_dimension", &args);
  sine_wave_ = *sine_wave;
  ASSERT(std::abs(sine_wave_.shift() - sine_wave_.default_shift()) < NEAR_ZERO,
    "use the intersection argument instead of shift in sine wave");
  sine_wave_.set_shift(intersection());
  solver_ = MakeSolverBrentDekker({{"tolerance", str(1e-8)}});
  // solver_ = MakeSolverBisection({{"tolerance", str(1e-8)}});
  // solver_ = MakeSolverNewtonRaphson({{"tolerance", str(1e-8)}, {"guess", "0"}});
  init_sine_dist_deriv_();
  check_all_used(args);
}

bool HalfSpaceSine::is_inside(const Position& point) const {
  const double e = point.coord(wave_dimension());
  const double f = point.coord(dimension());
  if (direction() == 1) {
    if (f >= sine_wave_.evaluate(e)) {
      return true;
    } else {
      return false;
    }
  } else {
    if (f >= sine_wave_.evaluate(e)) {
      return false;
    } else {
      return true;
    }
  }
}

double HalfSpaceSine::nearest_distance(const Position& point) const {
  TRACE("point " << point.str());
  TRACE("intersection: " << intersection());
  sine_dist_deriv_->set_point(point);
  const double e = point.coord(wave_dimension());
  TRACE(sine_wave().nearest_minimum(e));
  TRACE(sine_wave().nearest_maximum(e));
  const double max = sine_wave().nearest_maximum(e);
  const double min = sine_wave().nearest_minimum(e);
  solver_->set_lower(min);
  solver_->set_upper(max);
  solver_->set_guess(0.5*(max + min));
  const double x_min = solver_->root(sine_dist_deriv_.get());
//  solver_->set_guess(point.coord(dimension()));
//  const double x_min = solver_->root(&fx);
  TRACE("x_min " << x_min);
  const double y_min = sine_wave().evaluate(x_min);
  TRACE("y_min " << y_min);
  const double dx = point.coord(wave_dimension()) - x_min;
  const double dy = point.coord(dimension()) - y_min;
  double dist = std::sqrt(dx*dx + dy*dy);
  if (is_inside(point)) dist *= -1.;
  return dist;
}

bool HalfSpaceSine::is_inside(const Position& point,
    const double diameter) const {
  if (is_inside(point)) return true;
  const double dist = nearest_distance(point);
  TRACE("dist " << dist);
  if (dist > 0.5*diameter) return false;
  return true;
}

void HalfSpaceSine::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_half_space_(ostr);
  feasst_serialize_version(2345, ostr);
  feasst_serialize(wave_dimension_, ostr);
  feasst_serialize(solver_, ostr);
  feasst_serialize_fstobj(sine_wave_, ostr);
}

HalfSpaceSine::HalfSpaceSine(std::istream& istr) : HalfSpace(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2345 == version, version);
  feasst_deserialize(&wave_dimension_, istr);
  // HWH for unknown reasons, this doesn't work
  // feasst_deserialize(solver_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      solver_ = solver_->deserialize(istr);
    }
  }
  feasst_deserialize_fstobj(&sine_wave_, istr);
  init_sine_dist_deriv_();
}

}  // namespace feasst
