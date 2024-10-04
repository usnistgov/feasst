#include <cmath>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/solver_newton_raphson.h"
#include "math/include/formula.h"
#include "math/include/constants.h"

namespace feasst {

FEASST_MAPPER(SolverNewtonRaphson, argtype({{"guess", "0"}, {"tolerance", "0"}}));

SolverNewtonRaphson::SolverNewtonRaphson(argtype args) : Solver(args) {
  class_name_ = "SolverNewtonRaphson";
}

std::shared_ptr<Solver> SolverNewtonRaphson::create(std::istream& istr) const {
  return std::make_shared<SolverNewtonRaphson>(istr);
}

SolverNewtonRaphson::SolverNewtonRaphson(std::istream& istr)
  : Solver(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6937, "version mismatch: " << version);
}

void SolverNewtonRaphson::serialize(std::ostream& ostr) const {
  serialize_solver_(ostr);
  feasst_serialize_version(6937, ostr);
}

double SolverNewtonRaphson::root(Formula * formula) {
  double x = guess();
  double h = formula->evaluate(x)/formula->derivative(x);
  double h_prev = 0.;
  DEBUG(x);
  while (std::abs(h) >= tolerance()) {
    const double f = formula->evaluate(x);
    DEBUG("f " << f);
    const double fprime = formula->derivative(x);
    DEBUG("fprime " << fprime);
    h = f/fprime;
    DEBUG(h + h_prev);
    DEBUG("h " << h);
    x -= h;
    DEBUG(x);
    h_prev = h;
  }
  DEBUG("done");
  return x;
}

}  // namespace feasst
