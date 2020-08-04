#include <cmath>
#include "math/include/solver_newton_raphson.h"
#include "math/include/formula.h"
#include "math/include/constants.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

SolverNewtonRaphson::SolverNewtonRaphson(const argtype& args) : Solver(args) {
  class_name_ = "SolverNewtonRaphson";
}

class MapSolverNewtonRaphson {
 public:
  MapSolverNewtonRaphson() {
    auto solver = MakeSolverNewtonRaphson(
      {{"guess", "0"}, {"tolerance", "0"}});
    solver->deserialize_map()["SolverNewtonRaphson"] = solver;
  }
};

static MapSolverNewtonRaphson mapper_ = MapSolverNewtonRaphson();

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
