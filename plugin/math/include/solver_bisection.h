
#ifndef FEASST_MATH_SOLVER_BISECTION_H_
#define FEASST_MATH_SOLVER_BISECTION_H_

#include <map>
#include <string>
#include <memory>
#include <vector>
#include "math/include/solver.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Find roots for the continuous function f(x) within the interval [a, b]
  assuming that f(a)*f(b)<0,
 */
class SolverBisection : public Solver {
 public:
  explicit SolverBisection(argtype args = argtype());

  double root(Formula * formula) override;

  std::shared_ptr<Solver> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SolverBisection(std::istream& istr);
  virtual ~SolverBisection() {}

 private:
  // temporary
  double fa_, fb_;
};

inline std::shared_ptr<SolverBisection> MakeSolverBisection(
    const argtype &args = argtype()) {
  return std::make_shared<SolverBisection>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_SOLVER_BISECTION_H_
