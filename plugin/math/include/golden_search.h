
#ifndef FEASST_MATH_GOLDEN_SEARCH_H_
#define FEASST_MATH_GOLDEN_SEARCH_H_

#include <map>
#include <string>
#include <memory>
#include <vector>
#include "math/include/minimize.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 See Numerical Recipes 3rd Edition: The Art of Scientific Computing.
 */
class GoldenSearch : public Minimize {
 public:
  explicit GoldenSearch(argtype args = argtype());
  void bracket(double * lower, double * upper, Formula * formula) override;
  std::shared_ptr<Minimize> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit GoldenSearch(std::istream& istr);
  virtual ~GoldenSearch() {}
};

inline std::shared_ptr<GoldenSearch> MakeGoldenSearch(
    const argtype &args = argtype()) {
  return std::make_shared<GoldenSearch>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_GOLDEN_SEARCH_H_
