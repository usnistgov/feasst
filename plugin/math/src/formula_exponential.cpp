
#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/formula_exponential.h"

namespace feasst {

FEASST_MAPPER(FormulaExponential,);

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

FormulaExponential::FormulaExponential(argtype args)
  : FormulaExponential(&args) { feasst_check_all_used(args); }
FormulaExponential::FormulaExponential(argtype * args) : Formula(args) {
  class_name_ = "FormulaExponential";
  set_f0(dble("f0", args, 0.));
  set_A(dble("A", args, 1.));
  set_B(dble("B", args, 1.));
}

}  // namespace feasst
