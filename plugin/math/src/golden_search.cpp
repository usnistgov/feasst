#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/formula.h"
#include "math/include/golden_search.h"

namespace feasst {

GoldenSearch::GoldenSearch(argtype args) : Minimize(args) {
  class_name_ = "GoldenSearch";
}

class MapGoldenSearch {
 public:
  MapGoldenSearch() {
    auto solver = MakeGoldenSearch(
      {{"lower", "0"}, {"upper", "0"}, {"tolerance", "0"}});
    solver->deserialize_map()["GoldenSearch"] = solver;
  }
};

static MapGoldenSearch mapper_ = MapGoldenSearch();

std::shared_ptr<Minimize> GoldenSearch::create(std::istream& istr) const {
  return std::make_shared<GoldenSearch>(istr);
}

GoldenSearch::GoldenSearch(std::istream& istr)
  : Minimize(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3252, "version mismatch: " << version);
}

void GoldenSearch::serialize(std::ostream& ostr) const {
  serialize_solver_(ostr);
  feasst_serialize_version(3252, ostr);
}

// easy but not optimal number of function evaluations.
void GoldenSearch::bracket(double * a, double * b, Formula * formula) {
  const double tol = tolerance();
  *a = lower();
  *b = upper();
  feasst_sort(a, b);
  const double invphi = 2./(std::sqrt(5.) + 1.);
  DEBUG("invphi " << invphi);
  double c = *b - (*b - *a)*invphi;
  double d = *a + (*b - *a)*invphi;
  while (std::abs(c - d) > tol) {
    DEBUG("a: " << *a << " b: " << *b << " c: " << c << " d: " << d);
    if (formula->evaluate(c) < formula->evaluate(d)) {
      *b = d;
    } else {
      *a = c;
    }
    c = *b - (*b - *a)*invphi;
    d = *a + (*b - *a)*invphi;
  }
}

}  // namespace feasst
