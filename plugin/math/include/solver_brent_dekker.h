
#ifndef FEASST_MATH_SOLVER_BRENT_DEKKER_H_
#define FEASST_MATH_SOLVER_BRENT_DEKKER_H_

#include <string>
#include <memory>
#include <vector>
#include "math/include/solver.h"

namespace feasst {

/**
  https://en.wikipedia.org/wiki/Brent%27s_method
 */
class SolverBrentDekker : public Solver {
 public:
  SolverBrentDekker(const argtype& args = argtype());
  double root(Formula * formula) override;
  std::shared_ptr<Solver> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SolverBrentDekker(std::istream& istr);
  virtual ~SolverBrentDekker() {}
};

inline std::shared_ptr<SolverBrentDekker> MakeSolverBrentDekker(
    const argtype &args = argtype()) {
  return std::make_shared<SolverBrentDekker>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_SOLVER_BRENT_DEKKER_H_
