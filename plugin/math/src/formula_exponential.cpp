
#include "math/include/formula_exponential.h"
#include "utils/include/utils_io.h"
#include "utils/include/debug.h"

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
  feasst_deserialize_version(istr);
  feasst_deserialize(&f0_, istr);
  feasst_deserialize(&A_, istr);
  feasst_deserialize(&B_, istr);
}

void FormulaExponential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(1, ostr);
  serialize_formula_(ostr);
  feasst_serialize(f0_, ostr);
  feasst_serialize(A_, ostr);
  feasst_serialize(B_, ostr);
}

}  // namespace feasst
