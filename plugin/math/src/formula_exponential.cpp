
#include <cmath>
#include "math/include/formula_exponential.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapFormulaExponential {
 public:
  MapFormulaExponential() {
    FormulaExponential().deserialize_map()["FormulaExponential"] =
      std::make_shared<FormulaExponential>();
  }
};

static MapFormulaExponential mapper_ = MapFormulaExponential();

std::shared_ptr<Formula> FormulaExponential::create(std::istream& istr) const {
  return std::make_shared<FormulaExponential>(istr);
}

FormulaExponential::FormulaExponential(std::istream& istr)
  : Formula(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5837, "version mismatch: " << version);
  feasst_deserialize(&f0_, istr);
  feasst_deserialize(&A_, istr);
  feasst_deserialize(&B_, istr);
}

void FormulaExponential::serialize(std::ostream& ostr) const {
  serialize_formula_(ostr);
  feasst_serialize_version(5837, ostr);
  feasst_serialize(f0_, ostr);
  feasst_serialize(A_, ostr);
  feasst_serialize(B_, ostr);
}

double FormulaExponential::evaluate(const double x) const {
  return f0_*exp(A_*std::pow(x - x0(), B_));
}

FormulaExponential::FormulaExponential(const argtype& args) : Formula(args) {
  class_name_ = "FormulaExponential";
  args_.init(args);
  set_f0(args_.key("f0").dflt("0").dble());
  set_A(args_.key("A").dflt("1").dble());
  set_B(args_.key("B").dflt("1").dble());
}

}  // namespace feasst
