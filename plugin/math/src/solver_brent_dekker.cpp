#include <cmath>
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/constants.h"
#include "math/include/formula.h"
#include "math/include/utils_math.h"
#include "math/include/solver_brent_dekker.h"

namespace feasst {

FEASST_MAPPER(SolverBrentDekker,
              argtype({{"lower", "0"}, {"upper", "0"}, {"tolerance", "0"}}));

SolverBrentDekker::SolverBrentDekker(argtype args) : Solver(args) {
  class_name_ = "SolverBrentDekker";
}

std::shared_ptr<Solver> SolverBrentDekker::create(std::istream& istr) const {
  return std::make_shared<SolverBrentDekker>(istr);
}

SolverBrentDekker::SolverBrentDekker(std::istream& istr)
  : Solver(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1021, "version mismatch: " << version);
}

void SolverBrentDekker::serialize(std::ostream& ostr) const {
  serialize_solver_(ostr);
  feasst_serialize_version(1021, ostr);
}

double SolverBrentDekker::root(Formula * formula) {
  double a = lower();
  TRACE("a " << a);
  double b = upper();
  TRACE("b " << b);
  double fa = formula->evaluate(a);
  TRACE("fa " << fa);
  double fb = formula->evaluate(b);
  TRACE("fb " << fb);
  const double tol = tolerance();
  double absfa = std::abs(fa);
  double absfb = std::abs(fb);
  if (absfa < tol) return a;
  if (absfb < tol) return b;
  double s = -1, d = -1;
  ASSERT(fa*fb < 0, "fa: " << fa << " fb: " << fb << " fa*fb: " << fa*fb);
  if (absfa < absfb) {
    feasst_swap(&a, &b);
    feasst_swap(&fa, &fb);
    feasst_swap(&absfa, &absfb);
  }
  double c = a;
  double fc = fa;
  bool mflag = true;
  while (true) {
    if (absfb < tol) return b;
    if (absfa < tol) return a;
    if (std::abs(a - b) < tol) return 0.5*(a + b);
    if ( (fa != fc) && (fb != fc) ) {
      // inverse quadratic interpolation
      s = a*fb*fc/(fa - fb)/(fa - fc)
        + b*fa*fc/(fb - fa)/(fb - fc)
        + c*fa*fb/(fc - fa)/(fc - fb);
    } else {
      // secant
      s = b - fb*(b - a)/(fb - fa);
    }
    const double absbmc = std::abs(b - c);
    const double abssmb = std::abs(s - b);
    const double abscmd = std::abs(c - d);
    if ( (!is_in_interval(s, (3*a + b)/4, b)) ||
         (mflag && abssmb >= 0.5*absbmc) ||
         (!mflag && abssmb >= 0.5*abscmd) ||
         (mflag && absbmc < tol) ||
         (!mflag && abscmd < tol) ) {
      // bisection
      s = 0.5*(a + b);
      mflag = true;
    } else {
      mflag = false;
    }
    const double fs = formula->evaluate(s);
    const double absfs = std::abs(fs);
    d = c;
    c = b;
    fc = fb;
    if (fa*fs < 0) {
      b = s;
      fb = fs;
      absfb = absfs;
    } else {
      a = s;
      fa = fs;
      absfa = absfs;
    }
    if (absfa < absfb) {
      feasst_swap(&a, &b);
      feasst_swap(&fa, &fb);
      feasst_swap(&absfa, &absfb);
    }
  }
}

}  // namespace feasst
