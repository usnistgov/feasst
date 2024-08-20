
#ifndef FEASST_MATH_SOLVER_NEWTON_RAPHSON_H_
#define FEASST_MATH_SOLVER_NEWTON_RAPHSON_H_

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
class SolverNewtonRaphson : public Solver {
 public:
  explicit SolverNewtonRaphson(argtype args = argtype());
  double root(Formula * formula) override;
  std::shared_ptr<Solver> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SolverNewtonRaphson(std::istream& istr);
  virtual ~SolverNewtonRaphson() {}

 private:
};

inline std::shared_ptr<SolverNewtonRaphson> MakeSolverNewtonRaphson(
    const argtype &args = argtype()) {
  return std::make_shared<SolverNewtonRaphson>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_SOLVER_NEWTON_RAPHSON_H_
