#include <cmath>
#include "math/include/solver_bisection.h"
#include "math/include/formula.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

SolverBisection::SolverBisection(const argtype& args) : Solver(args) {
  class_name_ = "SolverBisection";
}

class MapSolverBisection {
 public:
  MapSolverBisection() {
    auto solver = MakeSolverBisection(
      {{"lower", "0"}, {"upper", "0"}, {"tolerance", "0"}});
    solver->deserialize_map()["SolverBisection"] = solver;
  }
};

static MapSolverBisection mapper_ = MapSolverBisection();

std::shared_ptr<Solver> SolverBisection::create(std::istream& istr) const {
  return std::make_shared<SolverBisection>(istr);
}

SolverBisection::SolverBisection(std::istream& istr)
  : Solver(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6937, "version mismatch: " << version);
}

void SolverBisection::serialize(std::ostream& ostr) const {
  serialize_solver_(ostr);
  feasst_serialize_version(6937, ostr);
}

double SolverBisection::root(Formula * formula) {
  double low = lower();
  double upp = upper();
  fa_ = formula->evaluate(low);
  fb_ = formula->evaluate(upp);
  if (std::abs(fa_) < tolerance()) return low;
  if (std::abs(fb_) < tolerance()) return upp;
  ASSERT(fa_*fb_ < 0, "fa: " << fa_ << " fb: " << fb_ << " fa*fb: " << fa_*fb_);
  while (upp - low >= tolerance()) {
    const double mid = 0.5*(upp + low);
    const double fc = formula->evaluate(mid);
    if (fc == 0) {
      return mid;
    } else if (fc*fa_ < 0) {
      upp = mid;
    } else {
      low = mid;
    }
  }
  return (upp + low)/2.;
}

}  // namespace feasst
